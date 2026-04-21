import gzip
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import copy_submission_file_to_gzip as script


class CopySubmissionFileToGzipTest(unittest.TestCase):
    def test_copies_plain_text_input_to_gzip_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "input.fa"
            output_path = tmpdir / "output.fa.gz"
            input_path.write_text(">seq1\nATGC\n")

            script.copy_to_gzip(input_path, output_path)

            with gzip.open(output_path, "rt") as handle:
                self.assertEqual(">seq1\nATGC\n", handle.read())

    def test_normalizes_gzipped_input_to_single_gzip_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "input.gff3.gz"
            output_path = tmpdir / "output.gff3.gz"
            with gzip.open(input_path, "wt") as handle:
                handle.write("##gff-version 3\n")

            script.copy_to_gzip(input_path, output_path)

            with gzip.open(output_path, "rt") as handle:
                self.assertEqual("##gff-version 3\n", handle.read())
