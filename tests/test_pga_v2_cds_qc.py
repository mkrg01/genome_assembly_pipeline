import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock
import json


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import pga_v2_cds_qc as script
from reference_cds_qc import CdsFeature
from organelle_annotation_utils import (
    extract_genbank_location_sequence,
    genbank_record_origin_sequence,
    reverse_complement_dna,
    split_genbank_records,
    validate_complete_cds_translation,
)

try:
    import Bio  # noqa: F401

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


CODON_PATTERN = [
    "GCT",
    "GAA",
    "TTC",
    "CCG",
    "AAG",
    "GGT",
    "CAA",
    "TGC",
    "ATT",
    "GAC",
]


def valid_cds(codon_count=80):
    middle = [
        CODON_PATTERN[index % len(CODON_PATTERN)]
        for index in range(codon_count - 2)
    ]
    return "ATG" + "".join(middle) + "TAA"


def origin_lines(sequence):
    lines = ["ORIGIN\n"]
    lower = sequence.lower()
    for offset in range(0, len(lower), 60):
        chunk = lower[offset : offset + 60]
        groups = " ".join(chunk[index : index + 10] for index in range(0, len(chunk), 10))
        lines.append(f"{offset + 1:9d} {groups}\n")
    lines.append("//\n")
    return lines


def genbank_text(record_id, sequence, feature_lines):
    return "".join(
        [
            f"LOCUS       {record_id:<16}{len(sequence):>7} bp    DNA     circular PLN 01-JAN-2000\n",
            "FEATURES             Location/Qualifiers\n",
            f"     source          1..{len(sequence)}\n",
            '                     /organism="Test species"\n',
            *feature_lines,
            *origin_lines(sequence),
        ]
    )


def cds_feature(location, gene, product="test protein", translation=None):
    lines = [
        f"     CDS             {location}\n",
        f'                     /gene="{gene}"\n',
        f'                     /product="{product}"\n',
        "                     /codon_start=1\n",
        "                     /transl_table=11\n",
    ]
    if translation is not None:
        lines.append(f'                     /translation="{translation}"\n')
    return lines


def exon_feature(location, gene):
    return [
        f"     exon            {location}\n",
        f'                     /gene="{gene}"\n',
    ]


def gene_feature(location, gene):
    return [
        f"     gene            {location}\n",
        f'                     /gene="{gene}"\n',
    ]


def write_reference(path, gene, cds_sequence):
    path.write_text(
        genbank_text(
            "ref",
            cds_sequence,
            cds_feature("1..240", gene),
        )
    )


def write_short_negative_target(path, cds_sequence):
    sequence = ["N"] * 900
    high = cds_sequence[:90]
    middle = cds_sequence[90:180]
    low = cds_sequence[180:240]
    placements = [
        ((600, 689), reverse_complement_dna(high)),
        ((300, 389), reverse_complement_dna(middle)),
        ((100, 159), reverse_complement_dna(low)),
    ]
    for (start, end), part in placements:
        self_length = end - start + 1
        if len(part) != self_length:
            raise AssertionError("fixture interval length mismatch")
        sequence[start - 1 : end] = list(part)
    full_sequence = "".join(sequence)
    path.write_text(
        genbank_text(
            "target",
            full_sequence,
            [
                *gene_feature("complement(100..689)", "clpP"),
                *exon_feature("complement(130..159)", "clpP"),
                *exon_feature("complement(300..389)", "clpP"),
                *exon_feature("complement(600..689)", "clpP"),
                *cds_feature(
                    "complement(join(130..159,300..389,600..689))",
                    "clpP",
                    translation="STALE",
                ),
            ],
        )
    )


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


class PgaV2FeatureBoundaryCoreTest(unittest.TestCase):
    def test_collinear_blocks_generate_negative_join_rescue_location(self):
        cds_sequence = valid_cds()
        sequence = ["N"] * 900
        for (start, end), part in [
            ((600, 689), reverse_complement_dna(cds_sequence[:90])),
            ((300, 389), reverse_complement_dna(cds_sequence[90:180])),
            ((100, 159), reverse_complement_dna(cds_sequence[180:240])),
        ]:
            sequence[start - 1 : end] = list(part)
        feature = CdsFeature(
            record_id="target",
            accession="",
            gene="clpP",
            product="",
            protein_id="",
            feature_index=4,
            location="complement(join(130..159,300..389,600..689))",
            strand="-",
            start=130,
            end=689,
            length=210,
            sequence="",
            codon_start=1,
            transl_table="11",
            protein="",
            protein_length=0,
            translation_source="",
            has_rna_editing_exception=False,
            validation_valid=False,
            validation_errors=[],
        )

        window = script.window_for_feature(feature, "".join(sequence), len(cds_sequence))
        blocks, reason = script.collinear_reference_blocks(
            cds_sequence,
            window["coding_sequence"],
        )
        intervals = [
            script.coding_block_to_genomic_interval(window, feature.strand, block)
            for block in blocks
        ]

        self.assertIsNone(reason)
        self.assertEqual(
            "complement(join(100..159,300..389,600..689))",
            script.cds_location_from_intervals(intervals, feature.strand),
        )


@unittest.skipUnless(BIOPYTHON_AVAILABLE, "Biopython is required for GenBank QC tests")
class PgaV2FeatureBoundaryRescueTest(unittest.TestCase):
    def test_feature_boundary_rescue_proposes_without_changing_annotation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reference_dir = root / "reference"
            reference_dir.mkdir()
            annotation = root / "target.gbk"
            cds_sequence = valid_cds()
            write_reference(reference_dir / "ref.gbk", "clpP", cds_sequence)
            write_short_negative_target(annotation, cds_sequence)
            before = annotation.read_text()

            summary = script.run_feature_boundary_rescue(
                annotation=annotation,
                reference_dir=reference_dir,
                warning_log=root / "missing_warning.log",
                qc_tsv=root / "boundary_qc.tsv",
                candidates_json=root / "boundary_candidates.json",
                fix_feature_boundaries=False,
            )

            self.assertFalse(summary["enabled"])
            self.assertEqual(1, summary["candidate_fix_count"])
            self.assertEqual(0, summary["applied_fix_count"])
            self.assertEqual(before, annotation.read_text())
            data = json.loads((root / "boundary_candidates.json").read_text())
            fix = data["candidates"][0]["candidate_fixes"][0]
            self.assertTrue(fix["auto_applicable"])
            self.assertEqual(
                "complement(join(100..159,300..389,600..689))",
                fix["new_location"],
            )

    def test_feature_boundary_rescue_applies_location_only_fix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reference_dir = root / "reference"
            reference_dir.mkdir()
            annotation = root / "target.gbk"
            cds_sequence = valid_cds()
            write_reference(reference_dir / "ref.gbk", "clpP", cds_sequence)
            write_short_negative_target(annotation, cds_sequence)
            original_origin = genbank_record_origin_sequence(
                split_genbank_records(annotation.read_text().splitlines(keepends=True))[0]
            )

            summary = script.run_feature_boundary_rescue(
                annotation=annotation,
                reference_dir=reference_dir,
                warning_log=root / "missing_warning.log",
                qc_tsv=root / "boundary_qc.tsv",
                candidates_json=root / "boundary_candidates.json",
                fix_feature_boundaries=True,
            )

            self.assertTrue(summary["enabled"])
            self.assertFalse(summary["sequence_modified"])
            self.assertEqual(1, summary["applied_fix_count"])
            self.assertEqual(1, summary["changed_cds_count"])
            self.assertEqual(3, summary["changed_exon_count"])
            text = annotation.read_text()
            self.assertIn(
                "CDS             complement(join(100..159,300..389,600..689))",
                text,
            )
            self.assertIn("exon            complement(100..159)", text)
            self.assertNotIn('/translation="STALE"', text)
            updated_origin = genbank_record_origin_sequence(
                split_genbank_records(text.splitlines(keepends=True))[0]
            )
            self.assertEqual(original_origin, updated_origin)
            rescued_cds = extract_genbank_location_sequence(
                "complement(join(100..159,300..389,600..689))",
                updated_origin,
            )
            self.assertEqual(cds_sequence, rescued_cds)
            self.assertTrue(
                validate_complete_cds_translation(rescued_cds, "11")["valid"]
            )

    def test_terminal_or_internal_codon_errors_only_are_not_candidates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reference_dir = root / "reference"
            reference_dir.mkdir()
            reference_cds = valid_cds()
            bad_start = "ACG" + reference_cds[3:]
            missing_stop = reference_cds[:-3] + "CAA"
            internal_stop = reference_cds[:30] + "TAA" + reference_cds[33:]
            reference_sequence = reference_cds * 3
            target_sequence = bad_start + missing_stop + internal_stop
            reference_features = [
                *cds_feature("1..240", "badStart"),
                *cds_feature("241..480", "missingStop"),
                *cds_feature("481..720", "internalStop"),
            ]
            target_features = [
                *cds_feature("1..240", "badStart"),
                *cds_feature("241..480", "missingStop"),
                *cds_feature("481..720", "internalStop"),
            ]
            (reference_dir / "ref.gbk").write_text(
                genbank_text("ref", reference_sequence, reference_features)
            )
            annotation = root / "target.gbk"
            annotation.write_text(genbank_text("target", target_sequence, target_features))

            summary = script.run_feature_boundary_rescue(
                annotation=annotation,
                reference_dir=reference_dir,
                warning_log=root / "missing_warning.log",
                qc_tsv=root / "boundary_qc.tsv",
                candidates_json=root / "boundary_candidates.json",
                fix_feature_boundaries=True,
            )

            self.assertEqual(0, summary["alignment_candidate_count"])
            self.assertEqual(0, summary["candidate_fix_count"])
            self.assertEqual(0, summary["applied_fix_count"])


if __name__ == "__main__":
    unittest.main()
