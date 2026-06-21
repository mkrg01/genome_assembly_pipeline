import sys
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from organelle_annotation_utils import format_post_curation_record


class PostCurationRecordTest(unittest.TestCase):
    def test_cds_auto_translation_findings_include_gene_and_location(self):
        text = format_post_curation_record(
            {
                "cds_auto_translation_validation": {
                    "cds_auto_translation_checked_cds_count": 12,
                    "cds_auto_translation_strict_table": "1",
                    "cds_auto_translation_comparison_table": "11",
                    "cds_auto_translation_strict_pass_count": 10,
                    "cds_auto_translation_table_compare_match_count": 12,
                    "cds_auto_translation_skipped_exception_count": 0,
                    "cds_auto_translation_skipped_partial_count": 0,
                    "cds_auto_translation_skipped_pseudo_count": 0,
                    "cds_auto_translation_blocking_issue_count": 2,
                    "cds_auto_translation_warning_count": 0,
                    "cds_auto_translation_passed": False,
                    "cds_auto_translation_issues": [
                        {
                            "severity": "error",
                            "record": "mt_ctg000001c",
                            "gene": "mttB",
                            "code": "strict_cds_translation_failed",
                            "message": (
                                "CDS is not valid for complete auto-translation "
                                "with table 1."
                            ),
                            "location": "186536..187324",
                        },
                        {
                            "severity": "error",
                            "record": "mt_ctg000001c",
                            "gene": "rpl16",
                            "code": "strict_cds_translation_failed",
                            "message": (
                                "CDS is not valid for complete auto-translation "
                                "with table 1."
                            ),
                            "location": "complement(159269..159706)",
                        },
                    ],
                }
            }
        )

        self.assertIn("CDS auto-translation QC: recorded 2", text)
        self.assertIn("mt_ctg000001c mttB (strict_cds_translation_failed)", text)
        self.assertIn("location=186536..187324", text)
        self.assertIn("mt_ctg000001c rpl16 (strict_cds_translation_failed)", text)
        self.assertIn("location=complement(159269..159706)", text)


if __name__ == "__main__":
    unittest.main()
