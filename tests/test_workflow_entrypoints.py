"""Check real Snakemake DAGs at the supported starting points."""

import gzip
import importlib.util
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


REPO = Path(__file__).resolve().parents[1]


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is required")
class WorkflowEntrypointsTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "config").mkdir()
        (self.root / "raw_data").mkdir()
        self.config = yaml.safe_load((REPO / "tests/workdir/config/config.yml").read_text())
        # Deliberately provide neither HiFi nor RNA-seq for downstream-only checks.
        self.config.update(external_assembly="raw_data/genome.fa", taxid=None)
        self.write("raw_data/genome.fa", ">chr1\n" + "ACGT" * 100 + "\n")
        self.write("raw_data/genes.gff3", "##gff-version 3\nchr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        self.write("raw_data/repeats.xm", "100 0 0 0 chr1 101 200 (200) + repeat#LTR/Gypsy 1 100 (0) 1\n")

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def run_workflow(self, target=None, entrypoint="Snakefile", dry_run=True):
        (self.root / "config/config.yml").write_text(yaml.safe_dump(self.config))
        cmd = [sys.executable, "-m", "snakemake", "--snakefile",
               str(REPO / "workflow" / entrypoint), "--directory", str(self.root),
               "--cores", "2", "--nolock"]
        if dry_run:
            cmd.append("--dry-run")
        if target:
            cmd.append(target)
        return subprocess.run(cmd, text=True, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=60)

    def jobs(self, target=None, entrypoint="Snakefile"):
        result = self.run_workflow(target, entrypoint)
        self.assertEqual(0, result.returncode, result.stdout)
        jobs = set(re.findall(r"^(?:local)?rule ([a-zA-Z0-9_]+):", result.stdout, re.M))
        # Compatibility files need a top-level default alias for included rules.
        return (jobs - {"default"}) | ({"all"} if "default" in jobs else set())

    def plot_inputs(self):
        self.config.update(external_gff3="raw_data/genes.gff3",
                           external_repeatmasker="raw_data/repeats.xm")

    def assert_no_assembly_jobs(self, jobs):
        self.assertFalse(jobs & {"hifiasm", "fastplong", "fcs_gx_screen", "oatk",
                                 "rename_assembly_by_length", "merqury"}, jobs)

    def test_external_annotation_matches_legacy_entrypoint_and_default(self):
        self.config["braker4_rnaseq_source"] = "varus"
        jobs = self.jobs("all")
        self.assertTrue({"stage_external_assembly", "repeatmasker", "braker4",
                         "format_for_submission", "plot_gene_synteny", "circos_plot",
                         "linear_plot", "busco_proteins_mode"} <= jobs, jobs)
        self.assert_no_assembly_jobs(jobs)
        self.assertFalse(any("organelle" in name for name in jobs), jobs)
        self.assertEqual(jobs, self.jobs())
        self.assertEqual(jobs, self.jobs("all", "Snakefile.annotation"))
        self.assertEqual(jobs, self.jobs(entrypoint="Snakefile.annotation"))
        self.assertEqual(jobs - {"all"}, self.jobs("annotation_all") - {"annotation_all"})

    def test_plots_with_external_annotations_skip_prediction_and_masking(self):
        self.plot_inputs()
        jobs = self.jobs("plots_all")
        self.assertTrue({"stage_external_assembly", "make_gene_bed", "make_repeat_bed",
                         "circos_plot", "linear_plot"} <= jobs, jobs)
        self.assert_no_assembly_jobs(jobs)
        self.assertFalse(jobs & {"braker4", "repeatmasker", "repeatmodeler", "omark",
                                 "download_busco_database", "prepare_gene_synteny"}, jobs)
        self.assertEqual(jobs - {"plots_all"},
                         self.jobs("all", "Snakefile.circos") - {"all"})
        self.assertEqual(self.jobs("all", "Snakefile.circos"),
                         self.jobs(entrypoint="Snakefile.circos"))

    def test_repeat_analysis_is_built_only_when_needed(self):
        self.plot_inputs()
        self.config["external_repeatmasker"] = None
        jobs = self.jobs("plots_all")
        self.assertIn("repeatmasker", jobs)
        self.assertNotIn("braker4", jobs)
        self.config["circos_plot_tracks"] = [self.config["circos_plot_tracks"][0]]
        jobs = self.jobs("plots_all")
        self.assertNotIn("repeatmasker", jobs)
        self.assertNotIn("tidk_search_for_circos_plot", jobs)

    def test_telomere_only_needs_no_annotation_or_repeat_analysis(self):
        self.config["circos_plot_tracks"] = [self.config["circos_plot_tracks"][-1]]
        jobs = self.jobs("plots_all")
        self.assertIn("tidk_search_for_circos_plot", jobs)
        self.assertFalse(jobs & {"braker4", "repeatmasker", "make_gene_bed"}, jobs)

    def test_missing_supplied_gff_fails_without_falling_back_to_braker(self):
        self.plot_inputs()
        self.config["external_gff3"] = "raw_data/missing.gff3"
        result = self.run_workflow("plots_all")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("raw_data/missing.gff3", result.stdout)
        self.assertIn("MissingInputException", result.stdout)

    def test_legacy_circos_accepts_pre_staged_pipeline_paths_without_reads(self):
        self.config["external_assembly"] = None
        for source, target in (
            ("raw_data/genome.fa", "results/renamed/assembly/primary/Testus_example.fa"),
            ("raw_data/genes.gff3", "results/braker4/primary/Testus_example/braker.gff3"),
            ("raw_data/repeats.xm", "results/repeatmasker/primary/Testus_example.fa.out.xm"),
        ):
            self.write(target, (self.root / source).read_text())
        jobs = self.jobs("all", "Snakefile.circos")
        self.assertIn("circos_plot", jobs)
        self.assertNotIn("braker4", jobs)
        self.assert_no_assembly_jobs(jobs)
        self.config["plots_only"] = True
        self.assertEqual(jobs - {"all"}, self.jobs("plots_all") - {"plots_all"})
        self.assertEqual(jobs, self.jobs("all"))
        self.assertEqual(jobs, self.jobs())

    def test_plot_only_mode_requires_gff_instead_of_running_gene_prediction(self):
        self.config["plots_only"] = True
        result = self.run_workflow("plots_all")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("MissingInputException", result.stdout)
        self.assertIn("braker.gff3", result.stdout)
        self.assertNotIn("No local RNA-Seq", result.stdout)

    def test_full_main_workflow_retains_assembly_and_organelle_paths(self):
        shutil.copytree(REPO / "tests/workdir/raw_data", self.root / "raw_data", dirs_exist_ok=True)
        self.config.update(external_assembly=None, taxid="122299")
        jobs = self.jobs("all")
        self.assertTrue({"hifiasm", "oatk", "fcs_gx_screen", "braker4",
                         "format_for_submission", "plot_gene_synteny"} <= jobs, jobs)
        self.assertEqual(jobs, self.jobs())

    def test_external_path_validation_and_annotation_requirement(self):
        for invalid in ("", True, ["raw_data/genome.fa"]):
            with self.subTest(invalid=invalid):
                self.config["external_assembly"] = invalid
                result = self.run_workflow("all")
                self.assertNotEqual(0, result.returncode)
                self.assertIn("non-empty path string", result.stdout)
        self.config["external_assembly"] = "raw_data/genome.fa"
        self.config["selected_assemblies"] = ["hap1", "hap2"]
        result = self.run_workflow("all")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("exactly one", result.stdout)
        self.config.update(external_assembly=None, selected_assemblies=["primary"])
        result = self.run_workflow("all", "Snakefile.annotation")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("'external_assembly' must be set", result.stdout)

    def test_stages_plain_and_gzipped_fasta_without_renaming_ids(self):
        for compressed in (False, True):
            with self.subTest(compressed=compressed):
                fasta = (self.root / "raw_data/genome.fa").read_bytes()
                self.config["external_assembly"] = "raw_data/genome.fa"
                if compressed:
                    (self.root / "raw_data/genome.fa.gz").write_bytes(gzip.compress(fasta))
                    self.config["external_assembly"] += ".gz"
                result = self.run_workflow("external_assembly_all", dry_run=False)
                self.assertEqual(0, result.returncode, result.stdout)
                self.assertEqual(fasta, (self.root / "results/external/assembly/primary/Testus_example.fa").read_bytes())

    def test_external_annotations_produce_bed_from_a_separate_workdir(self):
        self.plot_inputs()
        self.config["plots_only"] = True
        self.write("raw_data/gene models.gff3", (self.root / "raw_data/genes.gff3").read_text())
        self.config["external_gff3"] = "raw_data/gene models.gff3"
        for target, expected in (
            ("results/circos_plot/gene/primary/Testus_example_gene.bed",
             "chr1\t0\t100\tID=g1\t.\t+\n"),
            ("results/circos_plot/repeatmasker_repeat/primary/Testus_example_repeat_all.bed",
             "chr1\t100\t200\tLTR/Gypsy\t.\t+\n"),
        ):
            with self.subTest(target=target):
                result = self.run_workflow(target, dry_run=False)
                self.assertEqual(0, result.returncode, result.stdout)
                self.assertEqual(expected, (self.root / target).read_text())


if __name__ == "__main__":
    unittest.main()
