import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import rename_assembly_by_length as script


class RenameAssemblyByLengthTest(unittest.TestCase):
    def test_sorts_by_length_and_writes_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "input.fa"
            output_path = tmpdir / "renamed.fa"
            mapping_path = tmpdir / "mapping.tsv"
            input_path.write_text(
                ">short description one\n"
                "AA\n"
                ">long description two\n"
                "ACGT\n"
                "AC\n"
                ">medium\n"
                "NNNN\n"
            )

            script.rename_assembly(
                input_path,
                output_path,
                mapping_path,
                "scaffold",
                "longstitch",
            )

            self.assertEqual(
                ">scaffold1\nACGT\nAC\n"
                ">scaffold2\nNNNN\n"
                ">scaffold3\nAA\n",
                output_path.read_text(),
            )
            with mapping_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                ["long", "medium", "short"],
                [row["original_id"] for row in rows],
            )
            self.assertEqual(["6", "4", "2"], [row["length"] for row in rows])
            self.assertEqual(
                ["scaffold1", "scaffold2", "scaffold3"],
                [row["new_name"] for row in rows],
            )
            self.assertTrue(all(row["source_stage"] == "longstitch" for row in rows))

    def test_uses_original_id_to_break_length_ties(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "input.fa"
            output_path = tmpdir / "renamed.fa"
            mapping_path = tmpdir / "mapping.tsv"
            input_path.write_text(">zeta\nAAAA\n>alpha\nCCCC\n")

            script.rename_assembly(
                input_path,
                output_path,
                mapping_path,
                "contig",
                "fcs",
            )

            with mapping_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(["alpha", "zeta"], [row["original_id"] for row in rows])
            self.assertEqual(">contig1\nCCCC\n>contig2\nAAAA\n", output_path.read_text())

    def test_rejects_duplicate_record_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.fa"
            input_path.write_text(">duplicate first\nAA\n>duplicate second\nTT\n")

            with self.assertRaisesRegex(ValueError, "Duplicate FASTA record ID"):
                script.rename_assembly(
                    input_path,
                    Path(tmpdir) / "renamed.fa",
                    Path(tmpdir) / "mapping.tsv",
                    "contig",
                    "fcs",
                )

    def test_rejects_empty_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.fa"
            input_path.write_text("")

            with self.assertRaisesRegex(ValueError, "No FASTA records"):
                script.rename_assembly(
                    input_path,
                    Path(tmpdir) / "renamed.fa",
                    Path(tmpdir) / "mapping.tsv",
                    "contig",
                    "fcs",
                )


if __name__ == "__main__":
    unittest.main()
