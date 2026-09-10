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
sys.path.insert(0, str(REPO / "workflow/scripts"))
from purge_dups_support import fasta_lengths, filter_bed, read_cutoffs, resolve_cutoffs

if importlib.util.find_spec("snakemake"):
    from snakemake.io import expand
else:
    expand = None


def common(**overrides):
    config = {"pipeline_version": "test", "organism_name": "Testus_example", "genome_version": "test", "oatk_lineage": "magnoliopsida", **overrides}
    return runpy.run_path(str(REPO / "workflow/rules/common.smk"), init_globals={"config": config, "expand": expand})


class PurgeDupsConfigTest(unittest.TestCase):
    def test_routes_with_and_without_scaffolding(self):
        for enabled in (False, True):
            for longstitch in (False, True):
                for hic in (False, True):
                    with self.subTest(purge=enabled, longstitch=longstitch, hic=hic):
                        rules = common(
                            selected_assemblies=["primary", "hap1", "hap2"],
                            purge_dups_enabled=enabled, longstitch_enabled=longstitch,
                            hic_reads_r1=["r1.fq"] if hic else None,
                            hic_reads_r2=["r2.fq"] if hic else None,
                        )
                        for selected in ("primary", "hap1", "hap2"):
                            post = "purge_dups/assembly" if enabled else "fcs/assembly"
                            expected = f"results/{post}/{selected}/Testus_example.fa"
                            self.assertEqual(rules["post_fcs_assembly_path"]("Testus_example", selected), expected)
                            before_yahs = f"results/longstitch/assembly/{selected}/Testus_example.fa" if longstitch else expected
                            self.assertEqual(rules["pre_yahs_assembly_path"]("Testus_example", selected), before_yahs)
                            before_rename = f"results/yahs/assembly/{selected}/Testus_example.fa" if hic else before_yahs
                            self.assertEqual(rules["pre_rename_assembly_path"]("Testus_example", selected), before_rename)
                            if not hic and not longstitch:
                                self.assertEqual(rules["pre_rename_source_stage"](selected), post.split("/")[0])

    @unittest.skipUnless(expand, "Snakemake is required")
    def test_targets_follow_selected_assemblies(self):
        self.assertEqual(common()["purge_dups_targets"], [])
        for selected in (["primary"], ["hap2"], ["hap1", "hap2"], ["primary", "hap1", "hap2"]):
            with self.subTest(selected=selected):
                rules = common(selected_assemblies=selected, purge_dups_enabled=True)
                self.assertEqual(rules["purge_dups_targets"], selected)
                targets = rules["purge_dups_all_inputs"]("Testus_example")
                for assembly in selected:
                    for stage in ("fcs", "purge_dups"):
                        self.assertIn(f"results/{stage}/seqkit/{assembly}/Testus_example_seqkit_stats.tsv", targets)
                        self.assertIn(f"results/{stage}/length/{assembly}/Testus_example_length.pdf", targets)
                        self.assertIn(f"results/{stage}/gc_content/{assembly}/Testus_example_gc_content.pdf", targets)
                        self.assertIn(f"results/{stage}/busco_genome/{assembly}/BUSCO_Testus_example.fa", targets)
                        self.assertIn(f"results/{stage}/merqury/{assembly}/Testus_example.merqury.qv", targets)
                        self.assertIn(f"results/{stage}/dotplot/{assembly}/Testus_example_self_dotplot.pdf", targets)
                    self.assertIn(f"results/purge_dups/merqury/{assembly}/Testus_example.merqury.completeness.stats", targets)
                self.assertFalse(any("results/purge_dups/qc/" in path for path in targets))
                disabled = common(selected_assemblies=selected, purge_dups_enabled=False)
                self.assertEqual(disabled["purge_dups_targets"], [])
                self.assertFalse(any("results/purge_dups/" in path for path in disabled["purge_dups_all_inputs"]("Testus_example")))

    def test_invalid_configuration(self):
        for key, value in (
            ("purge_dups_enabled", "true"),
            ("purge_dups_cutoffs", [True, 20, 100]), ("purge_dups_cutoffs", [5, 100, 20]),
            ("purge_dups_cutoffs", "5,20,100"), ("purge_dups_cutoffs", [1, 2]),
        ):
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                common(**{key: value})


class PurgeDupsValidationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_fasta_and_bed_validation_preserves_depth_only_calls(self):
        fasta, bed, output = (self.root / name for name in ("assembly.fa", "raw.bed", "filtered.bed"))
        fasta.write_text(">primary\nACGTACGT\n>haplotig\nACGT\n")
        bed.write_text("primary\t0\t8\tHIGHCOV\nhaplotig\t0\t4\tHAPLOTIG\nprimary\t2\t6\tOVLP\n")
        filter_bed(fasta, bed, output)
        self.assertNotIn("HIGHCOV", output.read_text())
        self.assertIn("HAPLOTIG", output.read_text())
        self.assertIn("OVLP", output.read_text())
        bed.write_text("primary\t0\t9\tHAPLOTIG\n")
        with self.assertRaisesRegex(ValueError, "coordinates"):
            filter_bed(fasta, bed, output)
        for data in ("", ">x\n", ">x:1\nACGT\n", ">x\nAC\n>x\nGT\n"):
            fasta.write_text(data)
            with self.assertRaises(ValueError):
                fasta_lengths(fasta)

    def test_cutoffs_reject_failed_inference(self):
        self.assertEqual(read_cutoffs("5 19 19 20 20 100"), [5, 19, 19, 20, 20, 100])
        for value in ("", "0 0 0 0 0 0", "5 10 20 15 40 100", "not a histogram"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                read_cutoffs(value)


@unittest.skipUnless(shutil.which("get_seqs"), "purge_dups binaries are required")
class PurgeDupsExtractionTest(unittest.TestCase):
    def test_real_automatic_cutoff_inference(self):
        import math
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stat, output, metadata = (root / name for name in ("PB.stat", "cutoffs.txt", "cutoffs.json"))
            stat.write_text("".join(
                f"{depth}\t{int(100000 * math.exp(-((depth - 20) / 4) ** 2) + 200000 * math.exp(-((depth - 40) / 5) ** 2))}\n"
                for depth in range(501)
            ))
            resolve_cutoffs(stat, output, metadata)
            self.assertEqual(json.loads(metadata.read_text())["mode"], "auto")
            self.assertTrue(20 < read_cutoffs(output.read_text())[3] < 40)

    def test_real_get_seqs_keeps_internal_regions_and_short_remainders(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, bed = root / "input.fa", root / "dups.bed"
            fasta.write_text(">end\nAAAACCCCGGGG\n>internal\nAAAACCCCGGGG\n>hap\nACGT\n>short\nAC\n")
            bed.write_text("end\t0\t4\tOVLP\ninternal\t4\t8\tOVLP\nhap\t0\t4\tHAPLOTIG\n")
            subprocess.run(["get_seqs", "-e", "-c", "-l", "0", "-m", "0", "-g", "0", "-p", str(root / "out"), str(bed), str(fasta)], check=True, capture_output=True)
            self.assertEqual(fasta_lengths(root / "out.purged.fa"), {"end_1": 8, "internal_1": 12, "short_1": 2})
            self.assertEqual(sum(fasta_lengths(root / "out.hap.fa").values()), 8)


@unittest.skipUnless(importlib.util.find_spec("snakemake") and shutil.which("purge_dups") and shutil.which("minimap2"), "Snakemake and purge_dups/minimap2 are required")
class PurgeDupsWorkflowTest(unittest.TestCase):
    def test_real_mapping_purging_and_empty_bed_rerun(self):
        with tempfile.TemporaryDirectory(prefix="purge workflow ") as tmp:
            root = Path(tmp)
            (root / "workflow").symlink_to(REPO / "workflow", target_is_directory=True)
            fasta = root / "results/fcs/assembly/primary/Test.fa"
            reads = root / "results/hifi_reads/merged/Test_hifi_reads_curated.fastq.gz"
            fasta.parent.mkdir(parents=True)
            reads.parent.mkdir(parents=True)
            rng = random.Random(19)
            genome = "".join(rng.choices("ACGT", k=200000))
            fasta.write_text(">primary\n" + genome + "\n>haplotig\n" + genome[60000:100000] + "\n")
            import gzip
            with gzip.open(reads, "wt") as handle:
                for number, start in enumerate(range(0, len(genome) - 9999, 333)):
                    sequence = genome[start:start + 10000]
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            snakefile = root / "Snakefile"
            snakefile.write_text(
                "purge_dups_targets = ['primary']\npurge_dups_cutoffs = [2, 22, 100]\n"
                f"include: {str(REPO / 'workflow/rules/purge_dups.smk')!r}\n"
            )
            env = dict(os.environ, XDG_CACHE_HOME=str(root / "cache"))
            target = "results/purge_dups/assembly/primary/Test.fa"
            command = [sys.executable, "-m", "snakemake", "-s", str(snakefile), "--cores", "2", "--rerun-triggers", "mtime", "--", target]
            result = subprocess.run(command, cwd=root, env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(sum(fasta_lengths(root / target).values()), len(genome))
            self.assertIn("HAPLOTIG", (root / "results/purge_dups/bed/primary/Test.raw.bed").read_text())
            # A no-duplicate assembly must survive a repeat run; empty BED is valid.
            fasta.write_text(">primary\n" + genome + "\n")
            result = subprocess.run(command, cwd=root, env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((root / target).read_text(), fasta.read_text())
            self.assertEqual((root / "results/purge_dups/removed/primary/Test.fa").stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
