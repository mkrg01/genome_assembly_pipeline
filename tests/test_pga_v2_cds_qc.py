import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import pga_v2_cds_qc as script


class PgaV2CdsQcMappingTest(unittest.TestCase):
    def test_hifi_mapping_omits_unmapped_sam_records(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sam_path = Path(tmpdir) / "hifi_to_chloroplast.sam"

            with mock.patch.object(script.subprocess, "run") as run:
                cmd = script.map_hifi_reads(
                    Path("chloroplast.fa"),
                    Path("hifi.fastq.gz"),
                    sam_path,
                    threads=4,
                )

        self.assertIn("--sam-hit-only", cmd)
        self.assertEqual(
            [
                "minimap2",
                "-a",
                "--sam-hit-only",
                "-x",
                "map-hifi",
                "-t",
                "4",
                "chloroplast.fa",
                "hifi.fastq.gz",
            ],
            cmd,
        )
        run.assert_called_once()
        self.assertEqual(cmd, run.call_args.args[0])
        self.assertTrue(run.call_args.kwargs["check"])


if __name__ == "__main__":
    unittest.main()
