import difflib
import json
import re
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path

from organelle_annotation_utils import (
    reverse_complement_dna,
    validate_complete_cds_translation,
)
from reference_cds_qc import (
    SHORT_LENGTH_MISSING_BP,
    SHORT_LENGTH_RATIO,
    CdsFeature,
    parse_genbank_records_with_sequences,
    run_reference_cds_qc,
)


BASE_FLANK_BP = 300
SEED_LENGTHS = (30, 25, 20, 15)
READ_SUPPORT_MIN_SPANNING = 20
READ_SUPPORT_MIN_CORRECTED_FRACTION = 0.90
READ_SUPPORT_MAX_CURRENT_FRACTION = 0.05
READ_SUPPORT_MIN_MAPQ = 20
READ_SUPPORT_FLANK_BP = 20
READ_SUPPORT_INDEL_SLOP_BP = 8

@dataclass
class SequenceFix:
    record_id: str
    gene: str
    edit_type: str
    start: int
    end: int
    sequence: str
    reason: str
    read_support: dict = field(default_factory=dict)


def parse_pga_warnings(path: Path):
    warnings = {}
    if not path.exists():
        return warnings
    for line in path.read_text().splitlines():
        match = re.search(r"Warning:\s+(\S+)\s+\([^)]*\)\s+(.+)", line)
        if match:
            warnings.setdefault(match.group(1), []).append(match.group(2).strip())
    return warnings


def window_for_feature(feature: CdsFeature, record_sequence: str, reference_length: int):
    missing_bp = max(0, reference_length - feature.length)
    flank = BASE_FLANK_BP + missing_bp
    start = max(1, feature.start - flank)
    end = min(len(record_sequence), feature.end + flank)
    window = record_sequence[start - 1 : end]
    if feature.strand == "-":
        coding_window = reverse_complement_dna(window)
    else:
        coding_window = window
    return {
        "genomic_start": start,
        "genomic_end": end,
        "coding_sequence": coding_window,
    }


def best_seed_diagonal(reference_sequence: str, target_sequence: str):
    for seed_length in SEED_LENGTHS:
        diagonals = {}
        for ref_index in range(0, max(0, len(reference_sequence) - seed_length + 1), 3):
            seed = reference_sequence[ref_index : ref_index + seed_length]
            if "N" in seed:
                continue
            start = 0
            while True:
                target_index = target_sequence.find(seed, start)
                if target_index < 0:
                    break
                diagonal = target_index - ref_index
                diagonals[diagonal] = diagonals.get(diagonal, 0) + 1
                start = target_index + 1
        if diagonals:
            return max(diagonals.items(), key=lambda item: item[1])
    return None


def target_subsequence_for_reference(reference_sequence: str, target_sequence: str):
    diagonal = best_seed_diagonal(reference_sequence, target_sequence)
    if diagonal is None:
        return None
    start = max(0, diagonal[0])
    end = min(
        len(target_sequence),
        start + len(reference_sequence) + BASE_FLANK_BP,
    )
    return {
        "start": start,
        "end": end,
        "sequence": target_sequence[start:end],
        "seed_diagonal": diagonal[0],
        "seed_hit_count": diagonal[1],
    }


def difference_opcodes(reference_sequence: str, target_sequence: str):
    matcher = difflib.SequenceMatcher(
        None,
        reference_sequence,
        target_sequence,
        autojunk=False,
    )
    for tag, ref_start, ref_end, target_start, target_end in matcher.get_opcodes():
        if tag in {"insert", "delete", "replace"}:
            yield tag, ref_start, ref_end, target_start, target_end


def best_valid_deletion_span(
    reference_sequence,
    target_sequence,
    target_start,
    target_end,
    deletion_length,
    transl_table,
):
    first = max(0, target_start - 10)
    last = min(len(target_sequence) - deletion_length, target_end + 10)
    choices = []
    for deletion_start in range(first, last + 1):
        deletion_end = deletion_start + deletion_length
        validation = validate_corrected_cds(
            reference_sequence,
            target_sequence,
            "delete",
            deletion_start,
            deletion_end,
            "",
            transl_table,
        )
        if not validation["valid"]:
            continue
        deleted = target_sequence[deletion_start:deletion_end]
        homopolymer = len(set(deleted.upper())) == 1
        choices.append(
            (
                int(not homopolymer),
                abs(deletion_start - target_start),
                deletion_start,
                deletion_end,
                validation,
            )
        )
    if not choices:
        return None
    _score, _distance, deletion_start, deletion_end, validation = min(choices)
    return deletion_start, deletion_end, validation


def coding_index_to_genomic(window, feature: CdsFeature, coding_index: int):
    if feature.strand == "-":
        return window["genomic_end"] - coding_index
    return window["genomic_start"] + coding_index


def validate_corrected_cds(
    reference_sequence: str,
    target_subsequence: str,
    edit_type: str,
    edit_start: int,
    edit_end: int,
    edit_sequence: str,
    transl_table: str,
):
    if edit_type == "delete":
        corrected = target_subsequence[:edit_start] + target_subsequence[edit_end:]
    else:
        corrected = (
            target_subsequence[:edit_start]
            + edit_sequence
            + target_subsequence[edit_start:]
        )
    corrected_cds = corrected[: len(reference_sequence)]
    validation = validate_complete_cds_translation(corrected_cds, transl_table)
    length_delta_after = abs(len(corrected_cds) - len(reference_sequence))
    return {
        "valid": validation["valid"] and length_delta_after == 0,
        "errors": validation["errors"],
        "corrected_cds_length": len(corrected_cds),
        "corrected_cds_start": corrected_cds[:3],
        "corrected_cds_stop": corrected_cds[-3:] if len(corrected_cds) >= 3 else "",
    }


def find_frameshift_candidates(
    target_features,
    qc_rows,
    reference_choice_by_feature,
    target_records,
):
    records_by_id = {record["id"]: record for record in target_records}
    row_by_key = {
        (row.record_id, row.gene, row.location): row
        for row in qc_rows
    }
    candidates = []
    for feature in target_features:
        row = row_by_key.get((feature.record_id, feature.gene, feature.location))
        reference = reference_choice_by_feature.get(id(feature))
        if row is None or reference is None or not row.needs_alignment:
            continue
        record = records_by_id.get(feature.record_id)
        if record is None:
            continue
        window = window_for_feature(feature, record["sequence"], reference["length"])
        target_match = target_subsequence_for_reference(
            reference["sequence"],
            window["coding_sequence"],
        )
        candidate = {
            "record_id": feature.record_id,
            "gene": feature.gene,
            "target_location": feature.location,
            "target_length": feature.length,
            "reference_length": reference["length"],
            "reference_path": reference["path"],
            "window": {
                "genomic_start": window["genomic_start"],
                "genomic_end": window["genomic_end"],
            },
            "candidate_fixes": [],
            "skipped_reason": None,
        }
        if target_match is None:
            candidate["skipped_reason"] = "no_seed_match"
            candidates.append(candidate)
            continue
        candidate["target_match"] = {
            "coding_start": target_match["start"] + 1,
            "coding_end": target_match["end"],
            "seed_diagonal": target_match["seed_diagonal"],
            "seed_hit_count": target_match["seed_hit_count"],
        }
        for tag, ref_start, ref_end, target_start, target_end in difference_opcodes(
            reference["sequence"],
            target_match["sequence"],
        ):
            ref_length = ref_end - ref_start
            target_length = target_end - target_start
            if tag == "insert" or (tag == "replace" and target_length > ref_length):
                indel_length = target_length - ref_length
                if indel_length % 3 == 0:
                    continue
                deletion_span = best_valid_deletion_span(
                    reference["sequence"],
                    target_match["sequence"],
                    target_start,
                    target_end,
                    indel_length,
                    feature.transl_table,
                )
                if deletion_span is None:
                    continue
                deletion_start, deletion_end, validation = deletion_span
                coding_start = target_match["start"] + deletion_start
                coding_end = target_match["start"] + deletion_end - 1
                genomic_coords = [
                    coding_index_to_genomic(window, feature, index)
                    for index in range(coding_start, coding_end + 1)
                ]
                edit_start = min(genomic_coords)
                edit_end = max(genomic_coords)
                edit_sequence = record["sequence"][edit_start - 1 : edit_end]
                candidate["candidate_fixes"].append(
                    {
                        "edit_type": "delete",
                        "genomic_start": edit_start,
                        "genomic_end": edit_end,
                        "sequence": edit_sequence,
                        "coding_indel_sequence": target_match["sequence"][
                            deletion_start:deletion_end
                        ],
                        "indel_length": indel_length,
                        "validation": validation,
                        "reason": (
                            "target has a non-triplet insertion relative to "
                            "the reference CDS; deleting it restores a complete CDS"
                        ),
                    }
                )
            elif tag == "delete":
                indel_length = ref_length
                if indel_length % 3 == 0:
                    continue
                coding_position = target_match["start"] + target_start
                genomic_position = coding_index_to_genomic(
                    window,
                    feature,
                    coding_position,
                )
                missing_sequence = reference["sequence"][ref_start:ref_end]
                sequence_to_insert = (
                    reverse_complement_dna(missing_sequence)
                    if feature.strand == "-"
                    else missing_sequence
                )
                validation = validate_corrected_cds(
                    reference["sequence"],
                    target_match["sequence"],
                    "insert",
                    target_start,
                    target_end,
                    missing_sequence,
                    feature.transl_table,
                )
                if validation["valid"]:
                    candidate["candidate_fixes"].append(
                        {
                            "edit_type": "insert",
                            "genomic_start": genomic_position,
                            "genomic_end": genomic_position,
                            "sequence": sequence_to_insert,
                            "coding_indel_sequence": missing_sequence,
                            "indel_length": indel_length,
                            "validation": validation,
                            "reason": (
                                "target has a non-triplet deletion relative to "
                                "the reference CDS; inserting it restores a complete CDS"
                            ),
                        }
                    )
        candidates.append(candidate)
    return candidates


def cigar_ops(cigar):
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def reference_span(pos, cigar):
    ref_pos = pos
    end = pos - 1
    for length, op in cigar_ops(cigar):
        if op in {"M", "D", "N", "=", "X"}:
            end = ref_pos + length - 1
            ref_pos += length
    return pos, end


def classify_read_support(pos, cigar, candidate_fix):
    edit_type = candidate_fix["edit_type"]
    start = candidate_fix["genomic_start"]
    end = candidate_fix["genomic_end"]
    indel_length = candidate_fix["indel_length"]
    span_start, span_end = reference_span(pos, cigar)
    if span_start > start - READ_SUPPORT_FLANK_BP or span_end < end + READ_SUPPORT_FLANK_BP:
        return None

    corrected = False
    conflicting = False
    ref_pos = pos
    for length, op in cigar_ops(cigar):
        if op in {"M", "=", "X"}:
            ref_pos += length
            continue
        if op == "D":
            deletion_start = ref_pos
            deletion_end = ref_pos + length - 1
            overlaps = not (
                deletion_end < start - READ_SUPPORT_INDEL_SLOP_BP
                or deletion_start > end + READ_SUPPORT_INDEL_SLOP_BP
            )
            if overlaps:
                if edit_type == "delete" and abs(length - indel_length) <= 1:
                    corrected = True
                else:
                    conflicting = True
            ref_pos += length
            continue
        if op == "I":
            insertion_pos = ref_pos
            near = (
                start - READ_SUPPORT_INDEL_SLOP_BP
                <= insertion_pos
                <= end + READ_SUPPORT_INDEL_SLOP_BP
            )
            if near:
                if edit_type == "insert" and abs(length - indel_length) <= 1:
                    corrected = True
                else:
                    conflicting = True
            continue
        if op in {"N"}:
            ref_pos += length

    if corrected:
        return "corrected"
    if conflicting:
        return "ambiguous"
    return "current"


def map_hifi_reads(reference_fasta: Path, hifi_reads: Path, sam_path: Path, threads=1):
    sam_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "minimap2",
        "-a",
        "--sam-hit-only",
        "-x",
        "map-hifi",
        "-t",
        str(max(1, int(threads))),
        str(reference_fasta),
        str(hifi_reads),
    ]
    with sam_path.open("w") as handle:
        subprocess.run(cmd, check=True, stdout=handle)
    return cmd


def add_read_support(candidates, sam_path: Path):
    candidate_fixes = [
        fix
        for candidate in candidates
        for fix in candidate.get("candidate_fixes", [])
    ]
    for fix in candidate_fixes:
        fix["read_support"] = {
            "spanning_reads": 0,
            "corrected_reads": 0,
            "current_reads": 0,
            "ambiguous_reads": 0,
            "corrected_fraction": 0.0,
            "current_fraction": 0.0,
            "passes_auto_fix_thresholds": False,
        }
    if not candidate_fixes or not sam_path.exists():
        return

    with sam_path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            if flag & 0x4:
                continue
            mapq = int(fields[4])
            if mapq < READ_SUPPORT_MIN_MAPQ:
                continue
            pos = int(fields[3])
            cigar = fields[5]
            if cigar == "*":
                continue
            for fix in candidate_fixes:
                classification = classify_read_support(pos, cigar, fix)
                if classification is None:
                    continue
                support = fix["read_support"]
                support["spanning_reads"] += 1
                support[f"{classification}_reads"] += 1

    for fix in candidate_fixes:
        support = fix["read_support"]
        spanning = support["spanning_reads"]
        if spanning:
            support["corrected_fraction"] = support["corrected_reads"] / spanning
            support["current_fraction"] = support["current_reads"] / spanning
        support["passes_auto_fix_thresholds"] = (
            spanning >= READ_SUPPORT_MIN_SPANNING
            and support["corrected_fraction"] >= READ_SUPPORT_MIN_CORRECTED_FRACTION
            and support["current_fraction"] <= READ_SUPPORT_MAX_CURRENT_FRACTION
        )


def selected_sequence_fixes(candidates):
    fixes = []
    for candidate in candidates:
        passing = [
            fix
            for fix in candidate.get("candidate_fixes", [])
            if fix.get("read_support", {}).get("passes_auto_fix_thresholds")
        ]
        if len(passing) != 1:
            continue
        fix = passing[0]
        fixes.append(
            SequenceFix(
                record_id=candidate["record_id"],
                gene=candidate["gene"],
                edit_type=fix["edit_type"],
                start=fix["genomic_start"],
                end=fix["genomic_end"],
                sequence=fix["sequence"],
                reason=fix["reason"],
                read_support=fix.get("read_support", {}),
            )
        )
    return fixes


def apply_sequence_fixes(records, fixes):
    fixes_by_record = {}
    for fix in fixes:
        fixes_by_record.setdefault(fix.record_id, []).append(fix)
    updated_records = []
    for record in records:
        sequence = record["sequence"]
        record_fixes = sorted(
            fixes_by_record.get(record["id"], []),
            key=lambda fix: fix.start,
            reverse=True,
        )
        for fix in record_fixes:
            if fix.edit_type == "delete":
                observed = sequence[fix.start - 1 : fix.end]
                if observed.upper() != fix.sequence.upper():
                    raise RuntimeError(
                        "Refusing to apply PGA v2 sequence fix because the "
                        f"expected sequence at {record['id']}:{fix.start}.."
                        f"{fix.end} is {fix.sequence!r}, observed {observed!r}."
                    )
                sequence = sequence[: fix.start - 1] + sequence[fix.end :]
            elif fix.edit_type == "insert":
                sequence = sequence[: fix.start - 1] + fix.sequence + sequence[fix.start - 1 :]
            else:
                raise RuntimeError(f"Unsupported sequence fix type: {fix.edit_type}")
        updated = dict(record)
        updated["sequence"] = sequence
        updated_records.append(updated)
    return updated_records


def write_candidates_json(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def run_cds_qc(
    annotation: Path,
    reference_dir: Path,
    warning_log: Path,
    qc_tsv: Path,
    candidates_json: Path,
    hifi_reads: Path | None = None,
    hifi_reference_fasta: Path | None = None,
    fix_sequence_frameshifts=False,
    work_dir: Path | None = None,
    threads=1,
):
    target_records = parse_genbank_records_with_sequences(annotation)
    pga_warnings = parse_pga_warnings(warning_log)
    qc_result = run_reference_cds_qc(
        annotation=annotation,
        reference_dir=reference_dir,
        qc_tsv=qc_tsv,
        organelle="chloroplast",
        tool="pga_v2",
        phase="pre_rna_editing",
        pga_warnings=pga_warnings,
        default_transl_table="11",
    )
    target_features = qc_result["target_features"]
    qc_rows = qc_result["qc_rows"]
    reference_choice_by_feature = qc_result["reference_choice_by_feature"]
    candidates = find_frameshift_candidates(
        target_features,
        qc_rows,
        reference_choice_by_feature,
        target_records,
    )

    mapping_command = None
    if fix_sequence_frameshifts:
        if hifi_reads is None:
            raise RuntimeError(
                "organelle_reference_cds_qc.chloroplast.fix_hifi_frameshifts "
                "is true, but no HiFi reads were provided to run_pga_v2.py."
            )
        if hifi_reference_fasta is None:
            raise RuntimeError(
                "organelle_reference_cds_qc.chloroplast.fix_hifi_frameshifts "
                "is true, but no chloroplast FASTA was provided for HiFi "
                "read mapping."
            )
        sam_path = (work_dir or annotation.parent) / "pga_v2_cds_qc" / "hifi_to_chloroplast.sam"
        mapping_command = map_hifi_reads(
            hifi_reference_fasta,
            hifi_reads,
            sam_path,
            threads=threads,
        )
        add_read_support(candidates, sam_path)

    fixes = selected_sequence_fixes(candidates) if fix_sequence_frameshifts else []
    candidates_data = {
        "annotation": str(annotation),
        "reference_dir": str(reference_dir),
        "warning_log": str(warning_log),
        "thresholds": {
            "base_flank_bp": BASE_FLANK_BP,
            "short_length_ratio": SHORT_LENGTH_RATIO,
            "short_length_missing_bp": SHORT_LENGTH_MISSING_BP,
            "read_support_min_spanning": READ_SUPPORT_MIN_SPANNING,
            "read_support_min_corrected_fraction": READ_SUPPORT_MIN_CORRECTED_FRACTION,
            "read_support_max_current_fraction": READ_SUPPORT_MAX_CURRENT_FRACTION,
        },
        "fix_sequence_frameshifts": bool(fix_sequence_frameshifts),
        "hifi_reads": str(hifi_reads) if hifi_reads else None,
        "hifi_mapping_command": mapping_command,
        "qc_row_count": len(qc_rows),
        "alignment_candidate_count": len(candidates),
        "candidate_fix_count": sum(
            len(candidate.get("candidate_fixes", [])) for candidate in candidates
        ),
        "selected_fix_count": len(fixes),
        "selected_fixes": [asdict(fix) for fix in fixes],
        "candidates": candidates,
    }
    write_candidates_json(candidates_json, candidates_data)
    return {
        "pga_v2_cds_qc_tsv": str(qc_tsv),
        "pga_v2_cds_frameshift_candidates_json": str(candidates_json),
        "pga_v2_cds_qc_row_count": len(qc_rows),
        "pga_v2_cds_qc_alignment_candidate_count": len(candidates),
        "pga_v2_cds_qc_candidate_fix_count": candidates_data["candidate_fix_count"],
        "pga_v2_cds_qc_selected_fix_count": len(fixes),
        "pga_v2_cds_qc_selected_fixes": [asdict(fix) for fix in fixes],
    }
