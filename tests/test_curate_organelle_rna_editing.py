import sys
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import curate_organelle_rna_editing as script
from organelle_annotation_utils import GenBankFeatureBlock


class RnaEditingCdsNoteTest(unittest.TestCase):
    def test_rna_sequence_inference_omits_category_prefix(self):
        self.assertEqual(
            "similar to RNA sequence, mRNA (same species)",
            script.RNA_EDITING_INFERENCE,
        )

    def test_generic_rna_editing_does_not_add_translation_note(self):
        self.assertIsNone(
            script.rna_editing_cds_note(
                editing_effects=["nonsynonymous"],
                applied_decisions=["accept"],
            )
        )

    def test_start_gain_keeps_specific_note(self):
        self.assertEqual(
            script.RNA_EDITING_START_GAIN_NOTE,
            script.rna_editing_cds_note(
                editing_effects=["start_gain"],
                applied_decisions=["accept"],
            ),
        )

    def test_stop_gain_keeps_specific_note(self):
        self.assertEqual(
            script.RNA_EDITING_STOP_GAIN_NOTE,
            script.rna_editing_cds_note(
                editing_effects=["stop_gain"],
                applied_decisions=["accept"],
            ),
        )

    def test_reference_inferred_generic_editing_does_not_add_reference_note(self):
        self.assertIsNone(
            script.rna_editing_cds_note(
                editing_effects=["nonsynonymous"],
                applied_decisions=["reference_inferred"],
            )
        )

    def test_reference_inferred_start_and_stop_gain_keep_specific_notes(self):
        self.assertEqual(
            [
                script.RNA_EDITING_START_GAIN_NOTE,
                script.RNA_EDITING_STOP_GAIN_NOTE,
            ],
            script.rna_editing_cds_note(
                editing_effects=["start_gain", "stop_gain"],
                applied_decisions=["reference_inferred"],
            ),
        )

    def test_update_cds_block_omits_empty_rna_editing_note(self):
        block_lines = [
            "     CDS             1..9\n",
            '                     /gene="abc"\n',
            '                     /product="example protein"\n',
        ]

        updated = script.update_cds_block_lines(
            block_lines,
            translation="MKT",
            table_id="1",
            add_rna_editing_exception=True,
            rna_editing_note=script.rna_editing_cds_note(
                editing_effects=["nonsynonymous"],
                applied_decisions=["accept"],
            ),
            rna_editing_inferences=[script.RNA_EDITING_INFERENCE],
        )

        updated_text = "".join(updated)
        self.assertIn('/exception="RNA editing"', updated_text)
        self.assertIn("/translation=", updated_text)
        self.assertNotIn('/note="RNA editing is required for translation"', updated_text)

    def test_update_cds_block_replaces_undetermined_start_note(self):
        block_lines = [
            "     CDS             1..9\n",
            '                     /gene="abc"\n',
            '                     /note="start codon is not determined"\n',
            '                     /product="example protein"\n',
        ]

        updated = script.update_cds_block_lines(
            block_lines,
            translation="MKT",
            table_id="1",
            add_rna_editing_exception=True,
            rna_editing_note=script.RNA_EDITING_START_GAIN_NOTE,
            rna_editing_inferences=[script.RNA_EDITING_INFERENCE],
        )

        updated_text = "".join(updated)
        self.assertNotIn('/note="start codon is not determined"', updated_text)
        self.assertEqual(1, updated_text.count(f'/note="{script.RNA_EDITING_START_GAIN_NOTE}"'))

    def test_update_cds_block_deduplicates_existing_start_gain_note(self):
        block_lines = [
            "     CDS             1..9\n",
            '                     /gene="abc"\n',
            '                     /note="start codon is not determined"\n',
            f'                     /note="{script.RNA_EDITING_START_GAIN_NOTE}"\n',
            '                     /product="example protein"\n',
        ]

        updated = script.update_cds_block_lines(
            block_lines,
            translation="MKT",
            table_id="1",
            add_rna_editing_exception=True,
            rna_editing_note=script.RNA_EDITING_START_GAIN_NOTE,
            rna_editing_inferences=[script.RNA_EDITING_INFERENCE],
        )

        updated_text = "".join(updated)
        self.assertNotIn('/note="start codon is not determined"', updated_text)
        self.assertEqual(1, updated_text.count(f'/note="{script.RNA_EDITING_START_GAIN_NOTE}"'))

    def test_update_cds_block_replaces_undetermined_stop_note(self):
        block_lines = [
            "     CDS             1..9\n",
            '                     /gene="abc"\n',
            '                     /note="stop codon is not determined"\n',
            '                     /product="example protein"\n',
        ]

        updated = script.update_cds_block_lines(
            block_lines,
            translation="MKT",
            table_id="1",
            add_rna_editing_exception=True,
            rna_editing_note=script.RNA_EDITING_STOP_GAIN_NOTE,
            rna_editing_inferences=[script.RNA_EDITING_INFERENCE],
        )

        updated_text = "".join(updated)
        self.assertNotIn('/note="stop codon is not determined"', updated_text)
        self.assertEqual(1, updated_text.count(f'/note="{script.RNA_EDITING_STOP_GAIN_NOTE}"'))

    def test_update_cds_block_adds_unquoted_transl_table(self):
        block_lines = [
            "     CDS             1..9\n",
            '                     /gene="abc"\n',
            '                     /product="example protein"\n',
        ]

        updated = script.update_cds_block_lines(
            block_lines,
            translation="MKT",
            table_id="11",
            add_rna_editing_exception=True,
            rna_editing_note=None,
            rna_editing_inferences=[script.RNA_EDITING_INFERENCE],
        )

        updated_text = "".join(updated)
        self.assertIn("/transl_table=11", updated_text)
        self.assertNotIn('/transl_table="11"', updated_text)


class ReferenceInferenceTest(unittest.TestCase):
    def test_refseq_nucleotide_accession_uses_refseq_inference(self):
        self.assertEqual(
            "similar to AA sequence:RefSeq:NC_068654.1",
            script.reference_inference_value({"accession": "NC_068654.1"}),
        )

    def test_protein_id_is_preferred_over_record_accession_for_inference(self):
        self.assertEqual(
            "similar to AA sequence:RefSeq:YP_010567104.1",
            script.reference_inference_value(
                {
                    "accession": "NC_068654.1",
                    "protein_id": "YP_010567104.1",
                }
            ),
        )

    def test_protein_id_is_used_as_inference_fallback(self):
        self.assertEqual(
            "similar to AA sequence:RefSeq:YP_010567104.1",
            script.reference_inference_value({"protein_id": "YP_010567104.1"}),
        )

    def test_insd_accession_uses_insd_inference(self):
        self.assertEqual(
            "similar to AA sequence:INSD:OQ130029.1",
            script.reference_inference_value({"accession": "OQ130029.1"}),
        )

    def test_reference_site_evidence_note_uses_protein_id(self):
        self.assertEqual(
            "inferred from reference protein: RefSeq:YP_010567104.1",
            script.reference_site_evidence_note(
                {
                    "accession": "NC_068654.1",
                    "record_id": "NC_068654",
                    "protein_id": "YP_010567104.1",
                }
            ),
        )


class RnaEditingSiteFeatureTest(unittest.TestCase):
    def test_reference_inferred_site_uses_generic_and_reference_notes(self):
        site = {
            "genomic_pos": 42,
            "cds_pos": 1,
            "codon_index": 0,
            "codon_pos": 0,
            "genomic_base": "C",
            "edited_base": "T",
            "effect": "start_gain",
            "rna_depth": 0,
            "rna_edited_reads": 0,
            "rna_edit_fraction": 0,
            "dna_depth": 10,
            "dna_alt_fraction": 0,
        }

        accepted_site = script.make_accepted_site(
            site,
            decision="reference_inferred",
            reference={
                "accession": "NC_068654.1",
                "record_id": "NC_068654",
                "protein_id": "YP_010567104.1",
            },
        )
        blocks = script.new_rna_editing_site_feature_blocks(
            {"accepted_sites": [accepted_site]},
            existing_site_features=set(),
        )

        block_text = "".join(blocks[0])
        self.assertNotIn("closely related species", block_text)
        self.assertIn(f'/note="{script.RNA_EDITING_SITE_NOTE}"', block_text)
        self.assertIn(
            (
                '/note="inferred from reference protein: '
                'RefSeq:YP_010567104.1"'
            ),
            block_text,
        )

    def test_obsolete_reference_inferred_site_note_prevents_duplicate_feature(self):
        existing_block = GenBankFeatureBlock(
            key="misc_feature",
            start=0,
            end=2,
            lines=[
                "     misc_feature    42\n",
                (
                    '                     /note="C to U RNA editing inferred '
                    'from closely related species"\n'
                ),
            ],
            feature_index=1,
            location="42",
        )

        existing = script.existing_rna_editing_site_features([existing_block])
        blocks = script.new_rna_editing_site_feature_blocks(
            {
                "accepted_sites": [
                    {
                        "genomic_pos": 42,
                        "genomic_base": "C",
                    }
                ]
            },
            existing_site_features=existing,
        )

        self.assertEqual([], blocks)

    def test_reference_note_is_added_to_existing_site_feature(self):
        existing_block = GenBankFeatureBlock(
            key="misc_feature",
            start=0,
            end=2,
            lines=[
                "     misc_feature    42\n",
                f'                     /note="{script.RNA_EDITING_SITE_NOTE}"\n',
            ],
            feature_index=1,
            location="42",
        )
        existing, existing_lines = script.existing_rna_editing_site_feature_registry(
            [existing_block]
        )

        blocks = script.new_rna_editing_site_feature_blocks(
            {
                "accepted_sites": [
                    {
                        "genomic_pos": 42,
                        "genomic_base": "C",
                        "site_notes": [
                            script.RNA_EDITING_SITE_NOTE,
                            (
                                "inferred from reference protein: "
                                "RefSeq:YP_010567104.1"
                            ),
                        ],
                    }
                ]
            },
            existing_site_features=existing,
            existing_site_feature_lines=existing_lines,
        )

        block_text = "".join(existing_block.lines)
        self.assertEqual([], blocks)
        self.assertEqual(1, block_text.count(f'/note="{script.RNA_EDITING_SITE_NOTE}"'))
        self.assertEqual(
            1,
            block_text.count(
                (
                    '/note="inferred from reference protein: '
                    'RefSeq:YP_010567104.1"'
                )
            ),
        )

    def test_existing_reference_note_is_not_duplicated(self):
        existing_block = GenBankFeatureBlock(
            key="misc_feature",
            start=0,
            end=3,
            lines=[
                "     misc_feature    42\n",
                f'                     /note="{script.RNA_EDITING_SITE_NOTE}"\n',
                (
                    '                     /note="inferred from reference protein: '
                    'RefSeq:YP_010567104.1"\n'
                ),
            ],
            feature_index=1,
            location="42",
        )
        existing, existing_lines = script.existing_rna_editing_site_feature_registry(
            [existing_block]
        )

        script.new_rna_editing_site_feature_blocks(
            {
                "accepted_sites": [
                    {
                        "genomic_pos": 42,
                        "genomic_base": "C",
                        "site_notes": [
                            script.RNA_EDITING_SITE_NOTE,
                            (
                                "inferred from reference protein: "
                                "RefSeq:YP_010567104.1"
                            ),
                        ],
                    }
                ]
            },
            existing_site_features=existing,
            existing_site_feature_lines=existing_lines,
        )

        block_text = "".join(existing_block.lines)
        self.assertEqual(
            1,
            block_text.count(
                (
                    '/note="inferred from reference protein: '
                    'RefSeq:YP_010567104.1"'
                )
            ),
        )

    def test_obsolete_reference_site_note_is_replaced_on_existing_feature(self):
        existing_block = GenBankFeatureBlock(
            key="misc_feature",
            start=0,
            end=2,
            lines=[
                "     misc_feature    42\n",
                (
                    '                     /note="C to U RNA editing inferred '
                    'from closely related species"\n'
                ),
            ],
            feature_index=1,
            location="42",
        )
        existing, existing_lines = script.existing_rna_editing_site_feature_registry(
            [existing_block]
        )

        blocks = script.new_rna_editing_site_feature_blocks(
            {
                "accepted_sites": [
                    {
                        "genomic_pos": 42,
                        "genomic_base": "C",
                        "site_notes": [
                            script.RNA_EDITING_SITE_NOTE,
                            (
                                "inferred from reference protein: "
                                "RefSeq:YP_010567104.1"
                            ),
                        ],
                    }
                ]
            },
            existing_site_features=existing,
            existing_site_feature_lines=existing_lines,
        )

        block_text = "".join(existing_block.lines)
        self.assertEqual([], blocks)
        self.assertNotIn("closely related species", block_text)
        self.assertIn(f'/note="{script.RNA_EDITING_SITE_NOTE}"', block_text)
        self.assertIn(
            (
                '/note="inferred from reference protein: '
                'RefSeq:YP_010567104.1"'
            ),
            block_text,
        )

    def test_obsolete_reference_cds_evidence_note_is_replaced(self):
        existing_block = GenBankFeatureBlock(
            key="misc_feature",
            start=0,
            end=3,
            lines=[
                "     misc_feature    42\n",
                f'                     /note="{script.RNA_EDITING_SITE_NOTE}"\n',
                (
                    '                     /note="inferred from reference CDS '
                    'evidence: RefSeq:NC_068654.1"\n'
                ),
            ],
            feature_index=1,
            location="42",
        )
        existing, existing_lines = script.existing_rna_editing_site_feature_registry(
            [existing_block]
        )

        script.new_rna_editing_site_feature_blocks(
            {
                "accepted_sites": [
                    {
                        "genomic_pos": 42,
                        "genomic_base": "C",
                        "site_notes": [
                            script.RNA_EDITING_SITE_NOTE,
                            (
                                "inferred from reference protein: "
                                "RefSeq:YP_010567104.1"
                            ),
                        ],
                    }
                ]
            },
            existing_site_features=existing,
            existing_site_feature_lines=existing_lines,
        )

        block_text = "".join(existing_block.lines)
        self.assertNotIn("RefSeq:NC_068654.1", block_text)
        self.assertIn(
            (
                '/note="inferred from reference protein: '
                'RefSeq:YP_010567104.1"'
            ),
            block_text,
        )


if __name__ == "__main__":
    unittest.main()
