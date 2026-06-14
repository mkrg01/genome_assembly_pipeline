import gzip
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import resolve_oatk_minimum_kmer_coverage as script


class ResolveOatkMinimumKmerCoverageTest(unittest.TestCase):
    def test_accepts_manual_positive_integer(self):
        resolution = script.resolve_oatk_minimum_kmer_coverage("0250")

        self.assertEqual("manual", resolution["mode"])
        self.assertEqual(250, resolution["resolved_minimum_kmer_coverage"])
        self.assertIsNone(resolution["auto_multiplier"])

    def test_auto_uses_fixed_multiplier_and_nuclear_coverage(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            reads = tmpdir / "reads.fastq.gz"
            stats = tmpdir / "assembly_stats.tsv"

            with gzip.open(reads, "wt") as handle:
                handle.write("@read1\nAAAAAAAAAA\n+\nIIIIIIIIII\n")
                handle.write("@read2\nCCCCC\n+\nIIIII\n")
            stats.write_text(
                "file\tformat\tnum_seqs\tsum_len\tmin_len\n"
                "assembly.fa\tFASTA\t1\t10\t10\n"
            )

            resolution = script.resolve_oatk_minimum_kmer_coverage(
                "auto",
                hifi_reads=reads,
                nuclear_assembly_stats=stats,
            )

        self.assertEqual("auto", resolution["mode"])
        self.assertEqual(7, resolution["auto_multiplier"])
        self.assertEqual(15, resolution["hifi_read_bases"])
        self.assertEqual(10, resolution["nuclear_assembly_bases"])
        self.assertEqual(1.5, resolution["nuclear_coverage"])
        self.assertEqual(11, resolution["resolved_minimum_kmer_coverage"])

    def test_auto_requires_read_and_assembly_inputs(self):
        with self.assertRaisesRegex(ValueError, "requires both"):
            script.resolve_oatk_minimum_kmer_coverage("auto")

    def test_rejects_non_positive_manual_value(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            script.resolve_oatk_minimum_kmer_coverage("0")


if __name__ == "__main__":
    unittest.main()
