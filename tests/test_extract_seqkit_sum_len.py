import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import extract_seqkit_sum_len as script


class ExtractSeqkitSumLenTest(unittest.TestCase):
    def test_extracts_sum_len_by_header_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            stats_path = Path(tmpdir) / "stats.tsv"
            stats_path.write_text(
                "file\tformat\tnum_seqs\tsum_len\tmin_len\n"
                "assembly.fa\tFASTA\t42\t1,234,567\t100\n"
            )

            self.assertEqual(1234567, script.extract_sum_len(stats_path))

    def test_rejects_missing_sum_len_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            stats_path = Path(tmpdir) / "stats.tsv"
            stats_path.write_text("file\tformat\nassembly.fa\tFASTA\n")

            with self.assertRaisesRegex(ValueError, "sum_len"):
                script.extract_sum_len(stats_path)


if __name__ == "__main__":
    unittest.main()
