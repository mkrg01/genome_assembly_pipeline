import argparse
import configparser
import csv
import gzip
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow/scripts"))
import download_braker4
import filter_gff_by_fasta_ids
import rename_submission_gene_models
import run_braker4
import select_longest_cds


GFF3 = """##gff-version 3
chr1\tBRAKER4\tgene\t1\t100\t.\t+\t.\tID=geneA;
chr1\tBRAKER4\tmRNA\t1\t100\t.\t+\t.\tID=transcriptA;Parent=geneA;
chr1\tBRAKER4\tCDS\t1\t6\t.\t+\t0\tID=cdsA;Parent=transcriptA;
chr1\tBRAKER4\tmRNA\t1\t100\t.\t+\t.\tID=transcriptB;Parent=geneA;
chr1\tBRAKER4\tCDS\t1\t9\t.\t+\t0\tID=cdsB;Parent=transcriptB;
chr1\tBRAKER4\tgene\t200\t300\t.\t+\t.\tID=geneB;
chr1\tBRAKER4\tmRNA\t200\t300\t.\t+\t.\tID=independent;Parent=geneB;
chr1\tBRAKER4\tCDS\t200\t208\t.\t+\t0\tID=cdsC;Parent=independent;
"""
CDS = ">transcriptA\nATGAAA\n>transcriptB\nATGAAAAAA\n>independent\nATGAAAAAA\n"
AA = ">transcriptA\nMK\n>transcriptB\nMKK\n>independent\nMKK\n"
GTF = """chr1\tBRAKER4\tCDS\t1\t6\t.\t+\t0\tgene_id "geneA"; transcript_id "transcriptA";
chr1\tBRAKER4\tCDS\t1\t9\t.\t+\t0\tgene_id "geneA"; transcript_id "transcriptB";
chr1\tBRAKER4\tCDS\t200\t208\t.\t+\t0\tgene_id "geneB"; transcript_id "independent";
"""


class Braker4Test(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        source = self.root / "source"
        source.mkdir()
        (source / "Snakefile").write_text("rule all:\n    input: []\n")
        manifest = self.root / "source.json"
        manifest.write_text(json.dumps({"version": "v0.5.0-beta", "sha256": "test"}))
        for name in ("genome.fa", "proteins.fa"):
            (self.root / name).write_text(">chr1\nATGAAA\n")
        for lib in ("tissue.rep1", "tissue.rep2"):
            for mate in (1, 2):
                (self.root / f"{lib}_{mate}.fastq").write_text("@r\nATG\n+\nIII\n")
        lineage = self.root / "busco/lineages/embryophyta_odb12"
        lineage.mkdir(parents=True)
        (lineage / "dataset.cfg").write_text("name=embryophyta_odb12\n")
        self.args = argparse.Namespace(
            source=source, source_manifest=manifest,
            template=REPO / "workflow/config/braker4.ini",
            sample="Testus_example_primary", genome=self.root / "genome.fa",
            proteins=self.root / "proteins.fa", lineage="embryophyta_odb12",
            rnaseq_source="local", varus_genus="", varus_species="",
            rnaseq_r1=[self.root / f"tissue.rep{i}_1.fastq" for i in (1, 2)],
            rnaseq_r2=[self.root / f"tissue.rep{i}_2.fastq" for i in (1, 2)],
            busco_downloads=self.root / "busco", threads=2, mem_mb=32000,
            container_cache=self.root / "containers", output_dir=self.root / "results",
            dry_run=False,
        )

    def outputs(self, run):
        results = run / "output" / self.args.sample / "results"
        results.mkdir(parents=True)
        for name, contents in zip(run_braker4.CORE_FILES, (GTF, GFF3, CDS, AA)):
            with gzip.open(results / (name + ".gz"), "wt") as handle:
                handle.write(contents)
        (results / ".done").touch()
        return results

    def test_etp_configuration_and_unique_library_aliases(self):
        with mock.patch.dict(os.environ, {"BRAKER4_USE_COMPLEASM_HINTS": "1"}):
            run, _, cmd, env = run_braker4.prepare_run(self.args)
        with (run / "samples.csv").open() as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["genome"], row["genome_masked"])
        self.assertEqual(row["protein_fasta"], str(self.args.proteins))
        self.assertEqual("", row["sra_ids"])
        self.assertEqual("", row["bam_files"])
        r1 = row["fastq_r1"].split(":")
        r2 = row["fastq_r2"].split(":")
        self.assertEqual(2, len(r1))
        self.assertEqual(2, len({Path(p).name.split('.')[0] for p in r1}))
        for actual, expected in zip(r1 + r2, self.args.rnaseq_r1 + self.args.rnaseq_r2):
            self.assertEqual(Path(actual).resolve(), expected)
        ini = configparser.ConfigParser()
        ini.read(run / "config.ini")
        self.assertEqual(2, ini.getint("SLURM_ARGS", "cpus_per_task"))
        self.assertEqual(32000, ini.getint("SLURM_ARGS", "mem_of_node"))
        self.assertTrue(ini.getboolean("PARAMS", "skip_busco"))
        self.assertTrue(ini.getboolean("PARAMS", "no_cleanup"))
        self.assertFalse(ini.getboolean("PARAMS", "run_omark"))
        self.assertFalse(ini.getboolean("PARAMS", "use_varus"))
        self.assertEqual("", row["varus_genus"])
        self.assertEqual("", row["varus_species"])
        self.assertNotIn("BRAKER4_USE_COMPLEASM_HINTS", env)
        self.assertEqual("local", cmd[cmd.index("--executor") + 1])
        self.assertIn("mem_mb=32000", cmd)
        self.assertEqual("none", cmd[cmd.index("--workflow-profile") + 1])
        self.assertIn(str(self.root), cmd[cmd.index("--singularity-args") + 1])

    def use_varus(self):
        self.args.rnaseq_source = "varus"
        self.args.varus_genus = "Arabidopsis"
        self.args.varus_species = "thaliana"
        self.args.rnaseq_r1 = []
        self.args.rnaseq_r2 = []

    def test_varus_without_local_reads_and_query_changes(self):
        local, _, _, _ = run_braker4.prepare_run(self.args)
        self.use_varus()
        for path in self.root.glob("*.fastq"):
            path.unlink()
        run, identity, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(local, run)
        with (run / "samples.csv").open() as handle:
            row = next(csv.DictReader(handle))
        for field in ("fastq_r1", "fastq_r2", "sra_ids", "bam_files"):
            self.assertEqual("", row[field])
        self.assertEqual("Arabidopsis", row["varus_genus"])
        self.assertEqual("thaliana", row["varus_species"])
        self.assertEqual("varus", identity["rnaseq_source"])
        self.assertEqual([], identity["rnaseq_r1"])
        ini = configparser.ConfigParser()
        ini.read(run / "config.ini")
        self.assertTrue(ini.getboolean("PARAMS", "use_varus"))
        self.assertEqual("docker://katharinahoff/varus-notebook:v0.0.6", ini["containers"]["varus_image"])
        self.assertEqual(run, run_braker4.prepare_run(self.args)[0])
        self.args.varus_species = "lyrata"
        self.assertNotEqual(run, run_braker4.prepare_run(self.args)[0])

    def test_invalid_or_mixed_rnaseq_modes_are_rejected(self):
        self.args.varus_genus = "Arabidopsis"
        with self.assertRaisesRegex(ValueError, "cannot be supplied"):
            run_braker4.prepare_run(self.args)
        self.use_varus()
        self.args.rnaseq_r1 = [self.root / "tissue.rep1_1.fastq"]
        with self.assertRaisesRegex(ValueError, "do not also supply"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_r1 = []
        for invalid in ("", "thaliana;echo", "thaliana extra"):
            self.args.varus_species = invalid
            with self.assertRaisesRegex(ValueError, "valid genus and species"):
                run_braker4.prepare_run(self.args)
        self.args.rnaseq_source = "local"
        with self.assertRaisesRegex(ValueError, "matching, non-empty"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_source = "automatic"
        with self.assertRaisesRegex(ValueError, "must be 'local' or 'varus'"):
            run_braker4.prepare_run(self.args)

    def test_retries_resume_but_changed_evidence_gets_a_new_workdir(self):
        first, _, _, _ = run_braker4.prepare_run(self.args)
        checkpoint = first / "completed-training"
        checkpoint.touch()
        self.args.threads = 1
        second, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertEqual(first, second)
        self.assertTrue(checkpoint.exists())
        self.args.rnaseq_r1.pop()
        self.args.rnaseq_r2.pop()
        third, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(first, third)
        self.assertTrue(checkpoint.exists())

    def test_assemblies_have_separate_training_directories(self):
        first, _, _, _ = run_braker4.prepare_run(self.args)
        self.args.sample = "Testus_example_hap1"
        second, _, _, _ = run_braker4.prepare_run(self.args)
        self.assertNotEqual(first, second)

    def test_missing_pair_or_lineage_is_rejected(self):
        self.args.rnaseq_r2.pop()
        with self.assertRaisesRegex(ValueError, "matching"):
            run_braker4.prepare_run(self.args)
        self.args.rnaseq_r1.pop()
        self.args.lineage = "embryophyta_odb10"
        with self.assertRaisesRegex(ValueError, "odb12"):
            run_braker4.prepare_run(self.args)

    def test_exports_compressed_results_through_submission_processing(self):
        run, identity, cmd, _ = run_braker4.prepare_run(self.args)
        self.outputs(run)
        run_braker4.export_results(run, self.args.output_dir, self.args.sample, identity, cmd)
        exported = self.args.output_dir
        self.assertEqual({"genes": 2, "transcripts": 3}, json.loads((exported / "run.json").read_text())["counts"])
        representative = exported / "representative.fa"
        select_longest_cds.select_longest(exported / "braker.codingseq", exported / "braker.gff3", representative)
        self.assertEqual({"transcriptB", "independent"}, run_braker4.fasta_ids(representative))
        filtered = exported / "representative.gff3"
        with mock.patch.object(sys, "argv", ["filter", "--fasta", str(representative),
                                             "--gff3", str(exported / "braker.gff3"), "--output", str(filtered)]):
            filter_gff_by_fasta_ids.main()
        genes, transcripts = rename_submission_gene_models.build_identifier_maps(exported / "braker.gff3", "Tesex")
        rename_submission_gene_models.rewrite_gff3(filtered, exported / "submission.gff3", genes, transcripts)
        rename_submission_gene_models.rewrite_fasta(representative, exported / "submission.fa", transcripts)
        self.assertEqual({"Tesex_000001.t2", "Tesex_000002.t1"}, run_braker4.fasta_ids(exported / "submission.fa"))
        self.assertIn("Parent=Tesex_000001.t2", (exported / "submission.gff3").read_text())

    def test_bad_results_do_not_replace_previous_exports(self):
        run, identity, cmd, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        old = self.args.output_dir / "braker.gff3"
        old.write_text("previous annotation\n")
        with gzip.open(results / "braker.aa.gz", "wt") as handle:
            handle.write(">wrong_id\nMK\n")
        with self.assertRaisesRegex(ValueError, "FASTA IDs differ"):
            run_braker4.export_results(run, self.args.output_dir, self.args.sample, identity, cmd)
        self.assertEqual("previous annotation\n", old.read_text())
        self.assertFalse((self.args.output_dir / "run.json").exists())

    def test_missing_collected_file_forces_collection_again(self):
        run, _, _, _ = run_braker4.prepare_run(self.args)
        results = self.outputs(run)
        (results / "braker.aa.gz").unlink()
        resumed, _, command, _ = run_braker4.prepare_run(self.args)
        self.assertEqual(run, resumed)
        self.assertEqual("collect_results", command[command.index("--forcerun") + 1])

    def test_longest_cds_ties_use_input_order(self):
        cds, gff, output = (self.root / name for name in ("cds.fa", "gene.gff3", "longest.fa"))
        cds.write_text(CDS.replace("ATGAAA\n", "ATGAAAAAA\n"))
        gff.write_text(GFF3)
        select_longest_cds.select_longest(cds, gff, output)
        self.assertEqual({"transcriptA", "independent"}, run_braker4.fasta_ids(output))

    def test_archive_checksum_and_extraction(self):
        archive = self.root / "source.tar.gz"
        with tarfile.open(archive, "w:gz") as handle:
            data = b"rule all:\n    input: []\n"
            member = tarfile.TarInfo("BRAKER4-test/Snakefile")
            member.size = len(data)
            handle.addfile(member, io.BytesIO(data))
        destination = self.root / "installed"
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            download_braker4.install_archive(archive, destination, "wrong")
        self.assertFalse(destination.exists())
        digest = download_braker4.sha256_file(archive)
        download_braker4.install_archive(archive, destination, digest)
        self.assertTrue((destination / "Snakefile").exists())

    @unittest.skipUnless(os.environ.get("BRAKER4_TEST_SOURCE"), "Set BRAKER4_TEST_SOURCE to check the real upstream DAG")
    def test_pinned_upstream_etp_dry_run(self):
        self.args.source = Path(os.environ["BRAKER4_TEST_SOURCE"])
        self.args.dry_run = True
        _, _, command, env = run_braker4.prepare_run(self.args)
        result = subprocess.run(command, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(0, result.returncode, result.stdout)
        self.assertIn("run_genemark_etp", result.stdout)
        self.assertIn("convert_gtf_to_gff3", result.stdout)
        self.assertIn("collect_results", result.stdout)

    @unittest.skipUnless(os.environ.get("BRAKER4_TEST_SOURCE"), "Set BRAKER4_TEST_SOURCE to check the real upstream DAG")
    def test_pinned_upstream_varus_etp_dry_run(self):
        self.use_varus()
        self.args.source = Path(os.environ["BRAKER4_TEST_SOURCE"])
        self.args.dry_run = True
        # Exercise argument parsing too: the parent passes empty R1/R2 lists.
        argv = [str(REPO / "workflow/scripts/run_braker4.py")]
        for key, value in vars(self.args).items():
            argv.append("--" + key.replace("_", "-"))
            if isinstance(value, list):
                argv.extend(map(str, value))
            elif value is not True:
                argv.append(str(value))
        result = subprocess.run([sys.executable, *argv], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertEqual(0, result.returncode, result.stdout)
        self.assertIn("rule run_varus:", result.stdout)
        self.assertIn("rule run_genemark_etp:", result.stdout)
        self.assertIn("rule collect_results:", result.stdout)
        self.assertNotIn("rule hisat2_align:", result.stdout)
        self.assertNotIn("rule fastp:", result.stdout)


if __name__ == "__main__":
    unittest.main()
