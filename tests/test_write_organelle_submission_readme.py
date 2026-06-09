import argparse
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import write_organelle_submission_readme as script


class WriteOrganelleSubmissionReadmeTest(unittest.TestCase):
    def test_writes_mitochondrial_readme(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_readme = tmpdir / "README.md"
            args = argparse.Namespace(
                organelle="mitochondrion",
                annotation_tool="pmga",
                assembly_name="Testus_example",
                genome_version="v1.0",
                input_genome=Path("results/oatk/oatk/Testus_example.mito.ctg.fasta"),
                input_annotation=Path("results/organelle_annotation/mitochondrion/pmga/Testus_example/Testus_example.mitochondrion.annotation.gbk"),
                output_genome=Path("results/release/v1.0/organelle/mitochondrion/Testus_example_v1.0_mitochondrion.fa.gz"),
                output_annotation=Path("results/release/v1.0/organelle/mitochondrion/Testus_example_v1.0_mitochondrion.gbk.gz"),
                output_readme=output_readme,
            )

            script.write_readme(args)

            text = output_readme.read_text()
            self.assertIn("mitochondrial genome", text)
            self.assertIn("results/oatk/oatk/Testus_example.mito.ctg.fasta", text)
            self.assertIn("results/organelle_annotation/mitochondrion/pmga/Testus_example/Testus_example.mitochondrion.annotation.gbk", text)
            self.assertIn("generated with `pmga`", text)

    def test_writes_readme_without_annotation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_readme = tmpdir / "README.md"
            args = argparse.Namespace(
                organelle="chloroplast",
                annotation_tool=None,
                assembly_name="Testus_example",
                genome_version="v1.0",
                input_genome=Path("results/oatk/oatk/Testus_example.pltd.ctg.fasta"),
                input_annotation=None,
                output_genome=Path("results/release/v1.0/organelle/chloroplast/Testus_example_v1.0_chloroplast.fa.gz"),
                output_annotation=None,
                output_readme=output_readme,
            )

            script.write_readme(args)

            text = output_readme.read_text()
            self.assertIn("chloroplast genome", text)
            self.assertIn("results/oatk/oatk/Testus_example.pltd.ctg.fasta", text)
            self.assertIn("no annotation tool was configured", text)
            self.assertNotIn("_annotation.gbk.gz", text)
