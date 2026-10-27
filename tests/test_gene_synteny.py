import csv
import importlib.util
import json
import os
from pathlib import Path
import random
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'workflow/scripts'))
from gene_synteny import prepare, read_genes, filter_hits, normalize_collinearity, blocks, search


def fixture(root, copies=3, count=15):
    rng = random.Random(271)
    seqs = [''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=180)) for _ in range(count)]
    genome, aa, gff = root / 'genome.fa', root / 'protein.fa', root / 'genes.gff3'
    genome.write_text(''.join(f'>ctg{c}\n' + 'A' * 20000 + '\n' for c in range(copies)))
    with aa.open('w') as proteins, gff.open('w') as features:
        for c in range(copies):
            for i, seq in enumerate(seqs):
                # Third copy is inverted. Sequences are identical across loci.
                index = count - i - 1 if c == 2 else i
                start, end = 101 + index * 1000, 700 + index * 1000
                name = f'g{c}_{i}'
                strand = '-' if c == 2 else '+'
                features.write(f'ctg{c}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={name}\nctg{c}\t.\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={name}.t1;Parent={name}\n')
                proteins.write(f'>{name}.t1\n{seq}\n')
    return genome, aa, gff


def hit(a, b, score=100, ev='1e-30', start=1, end=180):
    return '\t'.join(map(str, [a, b, 100, end-start+1, 0, 0, start, end, start, end, ev, score, 180, 180])) + '\n'


class GenePreparationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = fixture(self.root)
        self.prepared = self.root / 'prepared'

    def test_preserves_identical_proteins_at_distinct_loci(self):
        prepare(*self.inputs, self.prepared)
        genes = read_genes(self.prepared / 'genes.tsv')
        self.assertEqual(len(genes), 45)
        self.assertEqual(len({g['gene_id'] for g in genes}), 45)
        self.assertEqual(genes[0]['start'], 100)
        self.assertEqual(genes[0]['end'], 700)
        self.assertEqual(genes[30]['gene_id'], 'g2_14')
        self.assertEqual(genes[30]['rank'], 0)
        self.assertIn('sc000000000\tsg000000001\t101\t700', (self.prepared / 'data.gff').read_text())

    def test_rejects_missing_id_or_parent_and_multiple_isoforms(self):
        original = self.inputs[2].read_text()
        for old, new in [('ID=g0_0.t1;', 'ID=missing;'), ('Parent=g0_0\n', 'Parent=absent\n')]:
            with self.subTest(new=new):
                self.inputs[2].write_text(original.replace(old, new))
                with self.assertRaises(ValueError):
                    prepare(*self.inputs, self.prepared)
        self.inputs[2].write_text(original + 'ctg0\t.\tmRNA\t101\t700\t.\t+\t.\tID=extra;Parent=g0_0\n')
        with self.inputs[1].open('a') as h:
            h.write('>extra\nMPEPTIDE\n')
        with self.assertRaisesRegex(ValueError, 'Multiple proteins'):
            prepare(*self.inputs, self.prepared)

    def test_rejects_coordinates_duplicates_and_invalid_proteins(self):
        original = self.inputs[2].read_text()
        self.inputs[2].write_text(original.replace('\t700\t', '\t20001\t', 1))
        with self.assertRaisesRegex(ValueError, 'coordinates'):
            prepare(*self.inputs, self.prepared)
        self.inputs[2].write_text(original)
        self.inputs[1].write_text('>g0_0.t1\nMP*EPTIDE\n')
        with self.assertRaisesRegex(ValueError, 'protein sequence'):
            prepare(*self.inputs, self.prepared)
        self.inputs[1].write_text('>g0_0.t1\nMP\n>g0_0.t1\nMP\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            prepare(*self.inputs, self.prepared)

    def test_compressed_and_percent_escaped_ids(self):
        import gzip
        self.inputs[0].write_text('>contig:one\n'+'A'*1000+'\n')
        self.inputs[1].write_text('>gene:1.t1\nMPEPTIDE*\n')
        text = 'contig%3Aone\t.\tgene\t1\t100\t.\t+\t.\tID=gene%3A1\ncontig%3Aone\t.\tmRNA\t1\t100\t.\t+\t.\tID=gene%3A1.t1;Parent=gene%3A1\n'
        path = self.root / 'genes.gff3.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
        prepare(self.inputs[0], self.inputs[1], path, self.prepared)
        self.assertEqual(read_genes(self.prepared/'genes.tsv')[0]['gene_id'], 'gene:1')


class HitAndBlockTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.raw, self.filtered = self.root/'raw', self.root/'filtered'
        self.genes = [dict(tool_id=n, gene_id=n, contig='chr'+n[0], start=i*100, end=i*100+50, rank=i, strand='+', protein_id=n) for i,n in enumerate(['a1','a2','b1','b2'])]

    def test_nonself_before_top_n_and_reciprocal_union(self):
        self.raw.write_text(hit('a1','a1',1000)+hit('a1','b1',100)+hit('a1','b2',50)+hit('b2','a1',99)+hit('b1','a1',90)+hit('a2','b2',100,end=20))
        stats = filter_hits(self.raw, self.genes, self.filtered, max_targets=1)
        text = self.filtered.read_text()
        self.assertEqual(stats['homology_pairs'],2)
        self.assertEqual(stats['self_hits'],1)
        self.assertNotIn('a1\ta1',text)
        self.assertIn('b2\ta1',text)
        self.assertNotIn('a2',text)

    def test_default_evalue_keeps_stringent_protein_filter(self):
        self.raw.write_text(hit('a1','b1',ev='1e-11')+hit('a2','b2',ev='1e-7'))
        stats=filter_hits(self.raw,self.genes,self.filtered)
        self.assertEqual(stats['evalue'],1e-10)
        self.assertEqual(stats['homology_pairs'],1)
        self.assertEqual(stats['failed_filters'],1)
        self.assertIn('a1\tb1',self.filtered.read_text())
        self.assertNotIn('a2\tb2',self.filtered.read_text())

    def test_invalid_hit_fails_instead_of_silently_disappearing(self):
        for data in ('bad\n',hit('unknown','a1'),hit('a1','b1',score=float('nan')),hit('a1','b1',start=0)):
            with self.subTest(data=data):
                self.raw.write_text(data)
                with self.assertRaises(ValueError):
                    filter_hits(self.raw,self.genes,self.filtered)

    def test_block_normalization_and_empty_output(self):
        path=self.root/'data.collinearity'
        path.write_text('## Alignment 0: score=100.0 e_value=1e-20 N=2 sc0&sc1 minus\n 0- 0: a1 b2 1e-30\n 0- 1: a2 b1 1e-20\n')
        block_rows, anchors=normalize_collinearity(path,self.genes)
        self.assertEqual(block_rows[0]['orientation'],'-')
        self.assertEqual(block_rows[0]['anchors'],2)
        self.assertEqual(anchors[0]['gene_a'],'a1')
        path.write_text('# No blocks\n')
        self.assertEqual(normalize_collinearity(path,self.genes),([],[]))
        path.write_text('## Alignment 0: score=100 e_value=1e-20 N=3 sc0&sc1 plus\n0- 0: a1 b1 1e-20\n')
        with self.assertRaisesRegex(ValueError,'Inconsistent'):
            normalize_collinearity(path,self.genes)


class ConfigurationTest(unittest.TestCase):
    def test_selected_assemblies_without_synteny_config(self):
        try:
            from snakemake.io import expand
        except ImportError:
            expand = None
        common=runpy.run_path(str(REPO/'workflow/rules/common.smk'),init_globals=dict(config=dict(pipeline_version='test',organism_name='Testus_example',genome_version='v1.0',oatk_lineage='magnoliopsida',selected_assemblies=['hap1','hap2']),expand=expand))
        if expand:
            paths=common['gene_synteny_all_inputs']('Test')
            self.assertEqual(len(paths),26)
            self.assertTrue(all('/primary/' not in p for p in paths))


@unittest.skipUnless(importlib.util.find_spec('matplotlib'), 'Matplotlib is required')
class PlotTest(unittest.TestCase):
    def test_modes_preserve_identical_gene_pairs_and_gene_free_contig(self):
        from plot_gene_synteny import coordinates
        genes=[dict(gene_id='a',contig='c',start=100,end=200,rank=0)]
        contigs=[dict(name='c',length=1000,genes=1),dict(name='empty',length=500,genes=0)]
        a, rows, total=coordinates(genes,contigs,'gene_order')
        b, _, total_bp=coordinates(genes,contigs,'mb')
        self.assertEqual(a,{'a':0.5})
        self.assertEqual(b,{'a':150})
        self.assertEqual(total,2)
        self.assertEqual(total_bp,1500)
        self.assertEqual(rows[-1][0],'empty')


@unittest.skipUnless(shutil.which('diamond') and shutil.which('MCScanX'), 'DIAMOND and MCScanX are required')
class RealToolTest(unittest.TestCase):
    def test_real_search_and_forward_reverse_blocks(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            prepare(*fixture(root),root/'prepared')
            raw=root/'diamond/hits.tsv'
            search(root/'prepared/proteins.fa',raw,2)
            blocks(root/'prepared',raw,root/'blocks')
            stats=json.loads((root/'blocks/summary.json').read_text())
            self.assertEqual(stats['evalue'],1e-10)
            diamond_command=stats['diamond']['command']
            self.assertEqual(float(diamond_command[diamond_command.index('--evalue')+1]),1e-10)
            self.assertEqual(stats['mcscanx']['command'], ['MCScanX', 'data', '-a'])
            self.assertEqual((stats['min_anchors'],stats['max_gaps'],stats['block_evalue']), (5,25,1e-5))
            self.assertGreaterEqual(stats['blocks'],3)
            self.assertGreaterEqual(stats['anchors'],45)
            with (root/'blocks/blocks.tsv').open() as h:
                rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual({r['orientation'] for r in rows},{'+','-'})
            self.assertGreaterEqual(stats['self_hits'],45)
            with (root/'blocks/anchors.bedpe').open() as handle:
                bedpe=list(csv.reader(handle,delimiter='\t'))
            self.assertEqual(len(bedpe),stats['anchors'])
            self.assertEqual(len(bedpe[0]),10)
            self.assertGreaterEqual(int(bedpe[0][1]),0)
            # Homology exists but fewer than five anchors: MCScanX itself returns no blocks.
            short=root/'short';short.mkdir()
            prepare(*fixture(short,copies=2,count=2),short/'prepared')
            search(short/'prepared/proteins.fa',short/'diamond/hits.tsv',2)
            blocks(short/'prepared',short/'diamond/hits.tsv',short/'blocks')
            short_stats=json.loads((short/'blocks/summary.json').read_text())
            self.assertGreater(short_stats['homology_pairs'],0)
            self.assertEqual(short_stats['blocks'],0)
            self.assertEqual((short/'blocks/anchors.bedpe').stat().st_size,0)
            # Real no-duplicate input: exact self hits alone must produce valid empty tables.
            solo=root/'solo';solo.mkdir()
            prepare(*fixture(solo,copies=1),solo/'prepared')
            search(solo/'prepared/proteins.fa',solo/'diamond/hits.tsv',2)
            blocks(solo/'prepared',solo/'diamond/hits.tsv',solo/'blocks')
            self.assertEqual(json.loads((solo/'blocks/summary.json').read_text())['blocks'],0)


@unittest.skipUnless(importlib.util.find_spec('snakemake') and importlib.util.find_spec('matplotlib') and shutil.which('diamond') and shutil.which('MCScanX'), 'Snakemake, Matplotlib, DIAMOND, and MCScanX are required')
class WorkflowTest(unittest.TestCase):
    def test_actual_four_stage_workflow_and_plot_only_rerun(self):
        with tempfile.TemporaryDirectory(prefix='gene synteny ') as temp:
            root=Path(temp)
            genome, protein, gff=fixture(root)
            for source, destination in [(genome,'results/repeatmasker/primary/Test.fa.masked'),(protein,'results/longest_cds/primary/Test_aa.fa'),(gff,'results/longest_cds/primary/Test.gff3')]:
                target=root/destination
                target.parent.mkdir(parents=True,exist_ok=True)
                target.symlink_to(source)
            snakefile=root/'Snakefile'
            snakefile.write_text('import os\nconfig = dict(organism_name="Testus_example", min_long_contig_length=10000)\n' + f'include: {str(REPO / "workflow/rules/gene_synteny.smk")!r}\n')
            # Rule scripts are resolved relative to the workflow entrypoint.
            (root/'scripts').symlink_to(REPO/'workflow/scripts',target_is_directory=True)
            (root/'envs').symlink_to(REPO/'workflow/envs',target_is_directory=True)
            targets=[f'results/gene_synteny/primary/Test/plots/{mode}.pdf' for mode in ('gene_order','mb')]
            command=[sys.executable,'-m','snakemake','-s',str(snakefile),'--cores','2','--rerun-triggers','mtime','--',*targets]
            result=subprocess.run(command,cwd=root,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            summary=json.loads((root/'results/gene_synteny/primary/Test/plots/plot.json').read_text())
            self.assertEqual(summary['figures'][0]['unique_anchor_pairs'],summary['figures'][1]['unique_anchor_pairs'])
            self.assertGreater(summary['figures'][0]['unique_anchor_pairs'],0)
            raw=root/'results/gene_synteny/primary/Test/diamond/hits.tsv'
            before=raw.stat().st_mtime_ns
            # Changing a drawing threshold must not rerun protein mapping.
            snakefile.write_text(snakefile.read_text().replace('min_long_contig_length=10000','min_long_contig_length=30000'))
            result=subprocess.run([sys.executable,'-m','snakemake','-s',str(snakefile),'--cores','2','--',*targets],cwd=root,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual(raw.stat().st_mtime_ns,before)
            summary=json.loads((root/'results/gene_synteny/primary/Test/plots/plot.json').read_text())
            self.assertEqual(summary['figures'][0]['contigs'],0)


if __name__=='__main__':
    unittest.main()
