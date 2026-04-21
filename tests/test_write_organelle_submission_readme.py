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
                organelle="mito",
                assembly_name="Testus_example",
                assembly_version="v1.0",
                input_genome=Path("results/oatk/oatk/Testus_example.mito.ctg.fasta"),
                input_annotation=Path("results/oatk/oatk/Testus_example.annot_mito.txt"),
                output_genome=Path("results/submission/organelle/mito/Testus_example_v1.0_mito_genome.fa.gz"),
                output_annotation=Path("results/submission/organelle/mito/Testus_example_v1.0_mito_annotation.txt.gz"),
                output_readme=output_readme,
            )

            script.write_readme(args)

            text = output_readme.read_text()
            self.assertIn("mitochondrial genome", text)
            self.assertIn("results/oatk/oatk/Testus_example.mito.ctg.fasta", text)
            self.assertIn("Oatk outputs selected by `oatk_organelle`", text)
