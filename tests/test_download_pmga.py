import io
import stat
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import download_pmga as script


def add_dir(tar, name, mode):
    info = tarfile.TarInfo(name)
    info.type = tarfile.DIRTYPE
    info.mode = mode
    tar.addfile(info)


def add_file(tar, name, content, mode):
    data = content.encode("utf-8")
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mode = mode
    tar.addfile(info, io.BytesIO(data))


class DownloadPmgaPermissionTest(unittest.TestCase):
    def test_unpack_archive_makes_extracted_bundle_deletable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            archive = tmpdir / "PMGA.tar.gz"
            bundle = tmpdir / "PMGA"

            with tarfile.open(archive, "w:gz") as tar:
                add_dir(tar, "PMGA", 0o555)
                add_dir(tar, "PMGA/usr", 0o555)
                add_dir(tar, "PMGA/usr/lib", 0o555)
                add_file(tar, "PMGA/usr/lib/tool.txt", "pmga", 0o444)

            script.unpack_archive(archive, bundle)

            self.assertTrue(bundle.exists())
            self.assertTrue((bundle / "usr" / "lib").stat().st_mode & stat.S_IWUSR)
            self.assertTrue((bundle / "usr" / "lib" / "tool.txt").stat().st_mode & stat.S_IWUSR)

            script.remove_tree(bundle)
            self.assertFalse(bundle.exists())

    def test_remove_tree_repairs_non_writable_directories_before_removing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bundle = Path(tmpdir) / "PMGA"
            nested = bundle / "usr" / "lib"
            nested.mkdir(parents=True)
            (nested / "tool.txt").write_text("pmga\n")
            nested.chmod(0o500)
            (bundle / "usr").chmod(0o500)
            bundle.chmod(0o500)

            script.remove_tree(bundle)

            self.assertFalse(bundle.exists())


if __name__ == "__main__":
    unittest.main()
