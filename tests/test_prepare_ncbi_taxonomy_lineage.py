import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import prepare_ncbi_taxonomy_lineage as script


class PrepareNcbiTaxonomyLineageTest(unittest.TestCase):
    def test_ensure_ete_cache_dir_creates_cache_under_home(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            home = Path(tmpdir)

            with mock.patch.object(script.Path, "home", return_value=home):
                cache_dir = script.ensure_ete_cache_dir()

            self.assertEqual(cache_dir, home / ".local" / "share" / "ete")
            self.assertTrue(cache_dir.is_dir())


if __name__ == "__main__":
    unittest.main()
