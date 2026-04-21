import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import rename_submission_gene_models as script


class RenameSubmissionGeneModelsTest(unittest.TestCase):
    def test_make_submission_prefix_requires_binomial_name(self):
        with self.assertRaisesRegex(ValueError, "at least genus and species"):
            script.make_submission_prefix("Arabidopsis")

    def test_rewrites_gff3_and_fasta_with_submission_ids(self):
        gff3_text = textwrap.dedent(
            """
            ##gff-version 3
            chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;
            chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1;
            chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=g1.t1.CDS1;Parent=g1.t1;
            chr1\tsrc\tgene\t200\t300\t.\t+\t.\tID=g2;
            chr1\tsrc\tmRNA\t200\t300\t.\t+\t.\tID=g2.t3;Parent=g2;
            """
        ).lstrip()
        fasta_text = textwrap.dedent(
            """
            >g1.t1 transcript one
            ATGC
            >g2.t3
            GCTA
            """
        ).lstrip()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_gff3 = tmpdir / "input.gff3"
            input_fasta = tmpdir / "input.fa"
            output_gff3 = tmpdir / "renamed.gff3"
            output_fasta = tmpdir / "renamed.fa"
            input_gff3.write_text(gff3_text)
            input_fasta.write_text(fasta_text)

            prefix = script.make_submission_prefix("Dioncophyllum_thollonii")
            gene_map, transcript_map = script.build_identifier_maps(input_gff3, prefix)
            feature_count = script.rewrite_gff3(
                input_gff3,
                output_gff3,
                gene_map,
                transcript_map,
            )
            record_count = script.rewrite_fasta(
                input_fasta,
                output_fasta,
                transcript_map,
            )

            self.assertEqual("Dioth", prefix)
            self.assertEqual(
                {"g1": "Dioth_000001", "g2": "Dioth_000002"},
                gene_map,
            )
            self.assertEqual(
                {
                    "g1.t1": "Dioth_000001.t1",
                    "g2.t3": "Dioth_000002.t3",
                },
                transcript_map,
            )
            self.assertEqual(5, feature_count)
            self.assertEqual(2, record_count)

            output_gff3_text = output_gff3.read_text()
            self.assertIn("ID=Dioth_000001;", output_gff3_text)
            self.assertIn("ID=Dioth_000001.t1;Parent=Dioth_000001;", output_gff3_text)
            self.assertIn(
                "ID=Dioth_000001.t1.CDS1;Parent=Dioth_000001.t1;",
                output_gff3_text,
            )
            self.assertIn("ID=Dioth_000002.t3;Parent=Dioth_000002;", output_gff3_text)

            output_fasta_text = output_fasta.read_text()
            self.assertIn(">Dioth_000001.t1 transcript one", output_fasta_text)
            self.assertIn(">Dioth_000002.t3", output_fasta_text)
