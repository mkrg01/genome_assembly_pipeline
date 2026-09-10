import hashlib
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "workflow/scripts"))
import dfam_md5_validator


class DfamChecksumTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.data = Path(self.tmp.name) / "database.h5.gz"
        self.data.write_bytes(b"test database" * 200000)
        self.checksum = self.data.with_suffix(".gz.md5")
        self.digest = hashlib.md5(self.data.read_bytes()).hexdigest()

    def test_published_checksum_formats(self):
        for suffix in ("", "  database.h5.gz", " *database.h5.gz"):
            with self.subTest(suffix=suffix):
                self.checksum.write_text(self.digest.upper() + suffix + "\n")
                dfam_md5_validator.check_md5(self.data, self.checksum)

    def test_corrupt_download_is_rejected(self):
        self.checksum.write_text(self.digest)
        self.data.write_bytes(b"truncated download")
        with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
            dfam_md5_validator.check_md5(self.data, self.checksum)

    def test_invalid_checksum_is_rejected(self):
        for value in ("", "<html>server error</html>", "a" * 31, "z" * 32):
            with self.subTest(value=value):
                self.checksum.write_text(value)
                with self.assertRaisesRegex(ValueError, "Invalid MD5"):
                    dfam_md5_validator.check_md5(self.data, self.checksum)

    def test_large_file_is_read_in_bounded_chunks(self):
        self.checksum.write_text(self.digest)

        class BoundedReader(io.BytesIO):
            def read(reader, size=-1):
                self.assertGreater(size, 0)
                self.assertLessEqual(size, 1024 * 1024)
                return super().read(size)

        reader = BoundedReader(self.data.read_bytes())
        original_open = open

        def open_file(path, *args, **kwargs):
            if Path(path) == self.data:
                return reader
            return original_open(path, *args, **kwargs)

        with mock.patch("builtins.open", side_effect=open_file):
            dfam_md5_validator.check_md5(self.data, self.checksum)


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is required")
class SoftmaskWorkflowTest(unittest.TestCase):
    """Exercise real rules with small, offline stand-ins for external tools."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="repeat tools ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        (self.root / "scripts").symlink_to(REPO / "workflow/scripts", target_is_directory=True)
        (self.root / "Testus_example.fa").write_text(">chr1\nACGTACGTACGT\n")
        self.snakefile = self.root / "Snakefile"
        self.snakefile.write_text(
            "import os, re\n"
            "config.setdefault('dfam_version', '4.0')\n"
            "config.setdefault('dfam_lineage_name', 'Viridiplantae')\n"
            "organism_name = 'Testus_example'\n"
            "selected_assembly_pattern = '(?:primary|hap1|hap2)'\n"
            "def downstream_assembly_path(name, selected): return name + '.fa'\n"
            f"include: {str(REPO / 'workflow/rules/softmask.smk')!r}\n"
        )
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ["PATH"],
                        XDG_CACHE_HOME=str(self.root / "cache"))
        self.tool("BuildDatabase", """
import sys
from pathlib import Path
prefix = sys.argv[sys.argv.index('-name') + 1]
for suffix in ('.00.nsq', '.nal', '.translation'):
    Path(prefix + suffix).write_text('split BLAST DB')
""")
        self.tool("RepeatModeler", """
import sys
from pathlib import Path
prefix = sys.argv[sys.argv.index('-database') + 1]
assert Path(prefix + '.00.nsq').is_file()
assert not Path(prefix + '.nsq').exists()
for suffix, data in (('-families.fa', '>model#DNA/hAT\\nACGT\\n'),
                     ('-families.stk', '# STOCKHOLM 1.0\\n//\\n'),
                     ('-rmod.log', 'completed\\n')):
    Path(prefix + suffix).write_text(data)
""")
        self.tool("wget", """
import gzip, hashlib, sys
from pathlib import Path
output = Path(sys.argv[sys.argv.index('-O') + 1])
url = sys.argv[-1]
assert '/Dfam_4.0/families/FamDB/dfam40.' in url
data = gzip.compress(b'fixture FamDB component', mtime=0)
output.write_bytes(hashlib.md5(data).hexdigest().encode() + b'  ' +
                   output.name.removesuffix('.md5').encode() if url.endswith('.md5') else data)
""")
        self.tool("python3", f"""
import os, sys
from pathlib import Path
if sys.argv[1] != '/opt/FamDB/famdb.py':
    os.execv({sys.executable!r}, [{sys.executable!r}, *sys.argv[1:]])
database = Path(sys.argv[sys.argv.index('-i') + 1])
assert sorted(p.name for p in database.glob('*.h5')) == [
    'dfam40.0.h5', 'dfam40.curated.consensus.0.h5',
    'dfam40.uncurated.consensus.0.h5', 'dfam40.uncurated.consensus.1.h5']
if 'info' in sys.argv:
    print('FamDB 3.0.0 / Dfam 4.0')
elif sys.argv[-1] != 'Missing taxon':
    print('>dfam#LTR/Gypsy\\nACGT')
""")
        self.tool("RepeatMasker", """
import sys
from pathlib import Path
assert int(sys.argv[sys.argv.index('-parallel') + 1]) >= 1
library = Path(sys.argv[sys.argv.index('-lib') + 1]).read_text()
assert '>model#DNA/hAT' in library and '>dfam#LTR/Gypsy' in library
assembly = Path(sys.argv[-1])
Path(assembly.name + '.masked').write_text(assembly.read_text().lower())
Path(assembly.name + '.out.xm').write_text('cross_match output\\n')
""")

    def tool(self, name, body):
        script = self.bin / name
        script.write_text(f"#!{sys.executable}\n" + body)
        script.chmod(0o755)

    def run_workflow(self, *args):
        return subprocess.run(
            [sys.executable, "-m", "snakemake", "--snakefile", str(self.snakefile),
             "--directory", str(self.root), "--cores", "1", "--nolock", *args],
            env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=60,
        )

    def test_consensus_download_and_split_database_through_masking(self):
        target = "results/repeatmasker/hap1/Testus_example.fa.masked"
        result = self.run_workflow(target)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("acgt", (self.root / target).read_text())
        families = self.root / "results/repeatmodeler/hap1/Testus_example-families.fa"
        self.assertTrue(families.is_file())
        database = self.root / "results/repeatmodeler/database/hap1/Testus_example"
        self.assertTrue((database / "Testus_example.00.nsq").is_file())
        self.assertFalse((database / "Testus_example-families.fa").exists())
        result = self.run_workflow(target)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("Nothing to be done", result.stdout)

    def test_empty_export_stops_the_workflow(self):
        result = self.run_workflow("export_dfam_repeat_fasta", "--config",
                                   "dfam_lineage_name=Missing taxon")
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("No Dfam consensus sequences exported",
                      (self.root / "logs/export_dfam_repeat_fasta.err").read_text())
        self.assertFalse((self.root / "results/repeatmasker/dfam/4.0/dfam_Missing taxon.repeat.fasta").exists())


if __name__ == "__main__":
    unittest.main()
