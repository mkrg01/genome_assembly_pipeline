import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import run_pmga as script

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    HAVE_BIO = True
except ImportError:
    HAVE_BIO = False


class PmgaExecutionTest(unittest.TestCase):
    def run_logged(self, root, text, code=0):
        def execute(command, *, stdout, stderr):
            stdout.write(text)
            self.assertEqual(stderr, subprocess.STDOUT)
            return subprocess.CompletedProcess(command, code)
        with patch.object(script.subprocess, 'run', side_effect=execute):
            return script.run_pmga('apptainer', 'image', root, root/'input.fa', '1', 'sample', root/'raw')

    def test_exit_zero_internal_exception_stops_and_retains_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for text in ['Traceback (most recent call last):\nValueError: bad\n',
                         'FileNotFoundError: missing.tbl\n',
                         'ValueError: Did not recognise the LOCUS line layout:\n']:
                with self.subTest(text=text):
                    with self.assertRaisesRegex(RuntimeError, 'despite exit code 0'):
                        self.run_logged(root, text)
                    self.assertEqual((root/'raw/pmga.execution.log').read_text(), text)

    def test_nonzero_exit_retains_diagnostic_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(RuntimeError, 'code 9'):
                self.run_logged(root, 'failed stage\n', 9)
            self.assertEqual((root/'raw/pmga.execution.log').read_text(), 'failed stage\n')

    def test_generic_error_word_does_not_fail_successful_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            command = self.run_logged(Path(tmp), 'error_rate = 0.01\ncompleted\n')
            self.assertIn('mgavas_m', command)


@unittest.skipUnless(HAVE_BIO, 'Biopython is required')
class PmgaOutputTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.raw = self.root/'raw'
        self.raw.mkdir()
        self.annotation = self.root/'canonical.gbk'
        self.annotation.write_text('previous successful annotation\n')
        self.record = {'id': 'contig', 'sequence': 'ATGAAATAA'}

    def write_raw(self, sequence='ATGAAATAA', identifier='sample'):
        r = SeqRecord(Seq(sequence), id=identifier, name=identifier, description='')
        r.annotations['molecule_type'] = 'DNA'
        SeqIO.write([r], self.raw/'sample.gbk', 'genbank')
        return self.raw/'sample.gbk'

    def test_source_only_output_is_allowed_when_format_id_and_bases_match(self):
        p = self.write_raw()
        self.assertEqual(script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample'), p)
        self.assertEqual(self.annotation.read_bytes(), p.read_bytes())

    def test_input_id_is_also_allowed(self):
        p = self.write_raw(identifier='contig')
        script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample')
        self.assertEqual(self.annotation.read_bytes(), p.read_bytes())

    def test_mismatch_and_malformed_locus_leave_existing_annotation_untouched(self):
        for sequence, identifier, malformed in [('ATGCCCTAA', 'sample', False),
                                                  ('ATGAAATAA', 'other', False),
                                                  ('ATGAAATAA', 'sample', True)]:
            with self.subTest(sequence=sequence, identifier=identifier, malformed=malformed):
                p = self.write_raw(sequence, identifier)
                if malformed:
                    lines = p.read_text().splitlines(keepends=True)
                    lines[0] = 'LOCUS       Triantha_japonica_mt_ctg000002c67041 bp   DNA\n'
                    p.write_text(''.join(lines))
                with self.assertRaisesRegex(RuntimeError, 'Invalid PMGA GenBank output'):
                    script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample')
                self.assertEqual(self.annotation.read_text(), 'previous successful annotation\n')

    def test_missing_output_fails_without_overwriting_annotation(self):
        with self.assertRaisesRegex(RuntimeError, 'No GenBank output'):
            script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample')
        self.assertEqual(self.annotation.read_text(), 'previous successful annotation\n')

    def test_multiple_records_are_not_silently_accepted(self):
        p = self.write_raw()
        p.write_text(p.read_text() * 2)
        with self.assertRaisesRegex(RuntimeError, 'Expected one GenBank record'):
            script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample')
        self.assertEqual(self.annotation.read_text(), 'previous successful annotation\n')

    def test_feature_locations_are_left_for_existing_post_curation(self):
        p = self.write_raw()
        text = p.read_text().replace('FEATURES             Location/Qualifiers\n',
                                    'FEATURES             Location/Qualifiers\n'
                                    '     gene            8..3\n'
                                    '                     /gene="abc"\n')
        p.write_text(text)
        script.copy_pmga_annotation(self.raw, self.annotation, self.record, 'sample')
        self.assertEqual(self.annotation.read_bytes(), p.read_bytes())
