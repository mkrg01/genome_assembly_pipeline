import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import out2xm


class OutToOutXmTest(unittest.TestCase):
    def test_parse_out_line_without_star(self):
        line = "123 4.5 0.0 1.2 contig1 10 50 (0) + RLG LTR/Gypsy 1 40 (0) 7"
        row = out2xm.parse_out_line_to_outxm(line)

        self.assertEqual(
            [
                "123",
                "4.5",
                "0.0",
                "1.2",
                "contig1",
                "10",
                "50",
                "(0)",
                "+",
                "RLG#LTR/Gypsy",
                "1",
                "40",
                "(0)",
            ],
            row,
        )

    def test_parse_out_line_with_star(self):
        line = "321 1.0 0.0 0.0 contig2 100 200 (30) C repeat1 DNA/hAT 3 90 (0) 8 *"
        row = out2xm.parse_out_line_to_outxm(line)

        self.assertEqual("*", row[-1])
        self.assertEqual("repeat1#DNA/hAT", row[9])

    def test_out_to_outxm_skips_headers_and_writes_rows(self):
        input_text = "\n".join(
            [
                "SW perc div. perc del. perc ins. query sequence",
                "123 4.5 0.0 1.2 contig1 10 50 (0) + RLG LTR/Gypsy 1 40 (0) 7",
                "321 1.0 0.0 0.0 contig2 100 200 (30) C repeat1 DNA/hAT 3 90 (0) 8 *",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_path = tmpdir / "repeatmasker.out"
            output_path = tmpdir / "repeatmasker.out.xm"
            input_path.write_text(input_text)

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                out2xm.out_to_outxm(str(input_path), str(output_path))

            output_lines = output_path.read_text().splitlines()
            self.assertEqual(2, len(output_lines))
            self.assertIn("RLG#LTR/Gypsy", output_lines[0])
            self.assertTrue(output_lines[1].endswith("\t*"))
            self.assertIn("wrote 2 rows", stderr.getvalue())
