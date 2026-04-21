import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import filter_gff_by_fasta_ids as script


class FilterGffByFastaIdsTest(unittest.TestCase):
    def test_filters_gff3_to_requested_transcripts(self):
        fasta_text = textwrap.dedent(
            """
            >g1.t1 representative transcript
            ATGC
            """
        ).lstrip()
        gff3_text = textwrap.dedent(
            """
            ##gff-version 3
            chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;
            chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1;
            chr1\tsrc\texon\t1\t40\t.\t+\t.\tID=g1.t1.exon1;Parent=g1.t1;
            chr1\tsrc\tgene\t200\t300\t.\t+\t.\tID=g2;
            chr1\tsrc\tmRNA\t200\t300\t.\t+\t.\tID=g2.t1;Parent=g2;
            chr1\tsrc\texon\t200\t260\t.\t+\t.\tID=g2.t1.exon1;Parent=g2.t1;
            """
        ).lstrip()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fasta_path = tmpdir / "keep.fa"
            gff3_path = tmpdir / "input.gff3"
            output_path = tmpdir / "filtered.gff3"
            fasta_path.write_text(fasta_text)
            gff3_path.write_text(gff3_text)

            keep_mrna_ids = script.parse_fasta_ids(fasta_path)
            keep_gene_ids, found_mrna_ids = script.collect_gene_ids(gff3_path, keep_mrna_ids)
            counts = script.write_filtered_gff3(
                gff3_path,
                output_path,
                found_mrna_ids,
                keep_gene_ids,
            )

            self.assertEqual({"g1.t1"}, keep_mrna_ids)
            self.assertEqual({"g1"}, keep_gene_ids)
            self.assertEqual({"g1.t1"}, found_mrna_ids)
            self.assertEqual((1, 1, 3), counts)

            output_text = output_path.read_text()
            self.assertIn("ID=g1;", output_text)
            self.assertIn("ID=g1.t1;Parent=g1;", output_text)
            self.assertIn("ID=g1.t1.exon1;Parent=g1.t1;", output_text)
            self.assertNotIn("ID=g2;", output_text)
            self.assertNotIn("ID=g2.t1;", output_text)

    def test_parse_fasta_ids_rejects_empty_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "empty.fa"
            fasta_path.write_text("")

            with self.assertRaisesRegex(ValueError, "No FASTA record IDs"):
                script.parse_fasta_ids(fasta_path)
