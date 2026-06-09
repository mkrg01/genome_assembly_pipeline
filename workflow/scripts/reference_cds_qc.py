#!/usr/bin/env python3
import argparse
import csv
import difflib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from organelle_annotation_utils import (
    extract_genbank_location_sequence,
    first_genbank_qualifier_value,
    genbank_record_origin_sequence,
    parse_genbank_feature_blocks,
    split_genbank_records,
    translate_complete_codons,
    validate_complete_cds_translation,
)


SHORT_LENGTH_RATIO = 0.90
SHORT_LENGTH_MISSING_BP = 60
REFERENCE_EXTENSIONS = {".gb", ".gbk", ".gbff"}


@dataclass
class CdsFeature:
    record_id: str
    accession: str
    gene: str
    product: str
    protein_id: str
    feature_index: int
    location: str
    strand: str
    start: int
    end: int
    length: int
    sequence: str
    codon_start: int
    transl_table: str
    protein: str
    protein_length: int
    translation_source: str
    validation_valid: bool
    validation_errors: list[str]


@dataclass
class CdsQcRow:
    organelle: str
    tool: str
    phase: str
    record_id: str
    gene: str
    product: str
    feature_index: int
    location: str
    strand: str
    query_cds_length: int
    query_protein_length: int
    query_translation_source: str
    reference_path: str
    reference_id: str
    reference_gene: str
    reference_product: str
    reference_cds_length: int | None
    reference_protein_length: int | None
    protein_identity: float | None
    protein_coverage: float | None
    query_coverage: float | None
    reference_coverage: float | None
    length_ratio: float | None
    missing_bp: int | None
    pga_warning: str
    needs_alignment: bool
    validation_valid: bool
    validation_errors: str
    review_reason: str
    manual_review_candidate: bool

    @property
    def target_length(self):
        return self.query_cds_length

    @property
    def reference_length(self):
        return self.reference_cds_length


QC_FIELDS = [
    "organelle",
    "tool",
    "phase",
    "record_id",
    "gene",
    "product",
    "feature_index",
    "location",
    "strand",
    "query_cds_length",
    "query_protein_length",
    "query_translation_source",
    "reference_path",
    "reference_id",
    "reference_gene",
    "reference_product",
    "reference_cds_length",
    "reference_protein_length",
    "protein_identity",
    "protein_coverage",
    "query_coverage",
    "reference_coverage",
    "length_ratio",
    "missing_bp",
    "pga_warning",
    "needs_alignment",
    "validation_valid",
    "validation_errors",
    "review_reason",
    "manual_review_candidate",
]

MANUAL_CANDIDATE_FIELDS = [
    "organelle",
    "tool",
    "phase",
    "record_id",
    "gene",
    "product",
    "feature_index",
    "location",
    "problem",
    "validation_errors",
    "reference_id",
    "reference_gene",
    "reference_product",
    "protein_identity",
    "protein_coverage",
    "length_ratio",
    "suggested_action",
]


def record_id_from_locus(record_lines):
    for line in record_lines:
        if line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 2:
                return parts[1]
    return "record"


def accession_from_record(record_lines):
    for line in record_lines:
        if line.startswith("VERSION"):
            parts = line.split()
            if len(parts) >= 2:
                return parts[1]
        if line.startswith("ACCESSION"):
            parts = line.split()
            if len(parts) >= 2:
                return parts[1]
    return ""


def normalized_product(product: str):
    return " ".join(product.lower().replace("-", " ").replace("_", " ").split())


def normalize_translation(protein: str):
    return "".join(protein.split()).rstrip("*")


def conceptual_translation(sequence: str, codon_start: int, transl_table: str):
    coding_sequence = sequence[codon_start - 1 :]
    protein, _ambiguous_codons = translate_complete_codons(
        coding_sequence,
        transl_table,
    )
    return normalize_translation(protein)


def parse_cds_features(path: Path, default_transl_table="1"):
    records = split_genbank_records(path.read_text().splitlines(keepends=True))
    features = []
    for record in records:
        record_id = record_id_from_locus(record)
        record_sequence = genbank_record_origin_sequence(record)
        _features_index, _origin_index, blocks = parse_genbank_feature_blocks(record)
        for block in blocks:
            if block.key != "CDS":
                continue
            genes = block.qualifiers.get("gene", [])
            if not genes or not block.intervals:
                continue
            location = block.location
            sequence = extract_genbank_location_sequence(location, record_sequence)
            codon_start_text = first_genbank_qualifier_value(
                block.lines,
                "codon_start",
                "1",
            )
            try:
                codon_start = int(codon_start_text)
            except ValueError:
                codon_start = 1
            if codon_start not in {1, 2, 3}:
                codon_start = 1
            transl_table = first_genbank_qualifier_value(
                block.lines,
                "transl_table",
                str(default_transl_table),
            )
            coding_sequence = sequence[codon_start - 1 :]
            validation = validate_complete_cds_translation(
                coding_sequence,
                transl_table,
            )
            translation = first_genbank_qualifier_value(
                block.lines,
                "translation",
                None,
            )
            if translation:
                protein = normalize_translation(translation)
                translation_source = "genbank_translation"
            else:
                protein = conceptual_translation(sequence, codon_start, transl_table)
                translation_source = "conceptual_translation"
            starts = [start for start, _end in block.intervals]
            ends = [end for _start, end in block.intervals]
            features.append(
                CdsFeature(
                    record_id=record_id,
                    accession=accession_from_record(record),
                    gene=genes[0],
                    product=first_genbank_qualifier_value(
                        block.lines,
                        "product",
                        "",
                    ),
                    protein_id=first_genbank_qualifier_value(
                        block.lines,
                        "protein_id",
                        "",
                    ),
                    feature_index=block.feature_index,
                    location=location,
                    strand=(
                        "-"
                        if location.replace(" ", "").startswith("complement(")
                        else "+"
                    ),
                    start=min(starts),
                    end=max(ends),
                    length=len(sequence),
                    sequence=sequence,
                    codon_start=codon_start,
                    transl_table=transl_table,
                    protein=protein,
                    protein_length=len(protein),
                    translation_source=translation_source,
                    validation_valid=validation["valid"],
                    validation_errors=validation["errors"],
                )
            )
    return features


def parse_genbank_records_with_sequences(path: Path):
    records = []
    for record in split_genbank_records(path.read_text().splitlines(keepends=True)):
        records.append(
            {
                "id": record_id_from_locus(record),
                "sequence": genbank_record_origin_sequence(record),
            }
        )
    return records


def reference_cds_by_gene(reference_dir: Path, default_transl_table="1"):
    by_gene = {}
    for path in sorted(reference_dir.iterdir()):
        if not path.is_file() or path.suffix.lower() not in REFERENCE_EXTENSIONS:
            continue
        for feature in parse_cds_features(path, default_transl_table):
            by_gene.setdefault(feature.gene, []).append(
                {
                    "path": str(path),
                    "record_id": feature.record_id,
                    "accession": feature.accession,
                    "gene": feature.gene,
                    "product": feature.product,
                    "protein_id": feature.protein_id,
                    "length": feature.length,
                    "sequence": feature.sequence,
                    "protein": feature.protein,
                    "protein_length": feature.protein_length,
                    "translation_source": feature.translation_source,
                    "location": feature.location,
                }
            )
    return by_gene


def reference_cds_by_product(references_by_gene):
    by_product = {}
    for references in references_by_gene.values():
        for reference in references:
            product = normalized_product(reference.get("product", ""))
            if product:
                by_product.setdefault(product, []).append(reference)
    return by_product


def protein_match_metrics(query_protein: str, reference_protein: str):
    if not query_protein or not reference_protein:
        return {
            "matches": 0,
            "protein_identity": None,
            "protein_coverage": None,
            "query_coverage": None,
            "reference_coverage": None,
        }
    matcher = difflib.SequenceMatcher(
        None,
        query_protein,
        reference_protein,
        autojunk=False,
    )
    matches = sum(block.size for block in matcher.get_matching_blocks())
    query_coverage = matches / len(query_protein)
    reference_coverage = matches / len(reference_protein)
    return {
        "matches": matches,
        "protein_identity": matches / max(len(query_protein), len(reference_protein)),
        "protein_coverage": min(query_coverage, reference_coverage),
        "query_coverage": query_coverage,
        "reference_coverage": reference_coverage,
    }


def reference_score(target: CdsFeature, reference):
    metrics = protein_match_metrics(target.protein, reference.get("protein", ""))
    identity = metrics["protein_identity"] or 0.0
    coverage = metrics["protein_coverage"] or 0.0
    length_delta = abs(reference["length"] - target.length)
    return (identity, coverage, -length_delta)


def closest_reference_feature(target: CdsFeature, references):
    if not references:
        return None
    return min(references, key=lambda item: abs(item["length"] - target.length))


def best_reference_feature(target: CdsFeature, references_by_gene, references_by_product):
    references = references_by_gene.get(target.gene, [])
    if not references:
        references = references_by_product.get(normalized_product(target.product), [])
    if not references:
        return None
    return max(references, key=lambda reference: reference_score(target, reference))


def needs_candidate_alignment(row: CdsQcRow):
    if row.reference_length is None or row.length_ratio is None or row.missing_bp is None:
        return False
    return (
        row.length_ratio < SHORT_LENGTH_RATIO
        or row.missing_bp >= SHORT_LENGTH_MISSING_BP
        or (bool(row.pga_warning) and row.missing_bp > 0)
    )


def row_review_reasons(row: CdsQcRow):
    reasons = []
    if not row.reference_id:
        reasons.append("no_reference_found")
    if not row.validation_valid:
        reasons.extend(row.validation_errors.split(";"))
    if row.length_ratio is not None and row.length_ratio < SHORT_LENGTH_RATIO:
        reasons.append("short_vs_reference")
    if row.missing_bp is not None and row.missing_bp >= SHORT_LENGTH_MISSING_BP:
        reasons.append("missing_bp_vs_reference")
    if row.protein_identity is not None and row.protein_identity < 0.70:
        reasons.append("low_protein_identity")
    if row.protein_coverage is not None and row.protein_coverage < 0.70:
        reasons.append("low_protein_coverage")
    if row.pga_warning:
        reasons.append("tool_warning")
    return reasons


def make_qc_rows(
    target_features,
    references_by_gene,
    pga_warnings=None,
    *,
    organelle,
    tool,
    phase,
):
    rows = []
    reference_choice_by_feature = {}
    references_by_product = reference_cds_by_product(references_by_gene)
    pga_warnings = pga_warnings or {}
    for feature in target_features:
        reference = best_reference_feature(
            feature,
            references_by_gene,
            references_by_product,
        )
        reference_choice_by_feature[id(feature)] = reference
        if reference is None:
            reference_length = None
            reference_protein_length = None
            length_ratio = None
            missing_bp = None
            metrics = {
                "protein_identity": None,
                "protein_coverage": None,
                "query_coverage": None,
                "reference_coverage": None,
            }
            reference_path = ""
            reference_id = ""
            reference_gene = ""
            reference_product = ""
        else:
            reference_length = reference["length"]
            reference_protein_length = reference.get("protein_length")
            length_ratio = feature.length / reference_length if reference_length else None
            missing_bp = reference_length - feature.length
            metrics = protein_match_metrics(
                feature.protein,
                reference.get("protein", ""),
            )
            reference_path = reference["path"]
            reference_id = reference["record_id"]
            reference_gene = reference["gene"]
            reference_product = reference.get("product", "")
        warning = "; ".join(pga_warnings.get(feature.gene, []))
        row = CdsQcRow(
            organelle=organelle,
            tool=tool,
            phase=phase,
            record_id=feature.record_id,
            gene=feature.gene,
            product=feature.product,
            feature_index=feature.feature_index,
            location=feature.location,
            strand=feature.strand,
            query_cds_length=feature.length,
            query_protein_length=feature.protein_length,
            query_translation_source=feature.translation_source,
            reference_path=reference_path,
            reference_id=reference_id,
            reference_gene=reference_gene,
            reference_product=reference_product,
            reference_cds_length=reference_length,
            reference_protein_length=reference_protein_length,
            protein_identity=metrics["protein_identity"],
            protein_coverage=metrics["protein_coverage"],
            query_coverage=metrics["query_coverage"],
            reference_coverage=metrics["reference_coverage"],
            length_ratio=length_ratio,
            missing_bp=missing_bp,
            pga_warning=warning,
            needs_alignment=False,
            validation_valid=feature.validation_valid,
            validation_errors=";".join(feature.validation_errors),
            review_reason="",
            manual_review_candidate=False,
        )
        row.needs_alignment = needs_candidate_alignment(row)
        reasons = row_review_reasons(row)
        row.review_reason = ";".join(reasons)
        row.manual_review_candidate = bool(reasons)
        rows.append(row)
    return rows, reference_choice_by_feature


def format_tsv_value(value):
    if isinstance(value, float):
        return f"{value:.6f}"
    if value is None:
        return ""
    return value


def write_qc_tsv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=QC_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: format_tsv_value(getattr(row, field))
                    for field in QC_FIELDS
                }
            )


def manual_candidate_rows(rows):
    candidates = []
    for row in rows:
        if row.validation_valid:
            continue
        candidates.append(
            {
                "organelle": row.organelle,
                "tool": row.tool,
                "phase": row.phase,
                "record_id": row.record_id,
                "gene": row.gene,
                "product": row.product,
                "feature_index": row.feature_index,
                "location": row.location,
                "problem": row.validation_errors,
                "validation_errors": row.validation_errors,
                "reference_id": row.reference_id,
                "reference_gene": row.reference_gene,
                "reference_product": row.reference_product,
                "protein_identity": format_tsv_value(row.protein_identity),
                "protein_coverage": format_tsv_value(row.protein_coverage),
                "length_ratio": format_tsv_value(row.length_ratio),
                "suggested_action": "manual_review",
            }
        )
    return candidates


def write_manual_candidates_tsv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANUAL_CANDIDATE_FIELDS,
            delimiter="\t",
        )
        writer.writeheader()
        for row in manual_candidate_rows(rows):
            writer.writerow(row)


def write_qc_json(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def run_reference_cds_qc(
    *,
    annotation: Path,
    reference_dir: Path,
    qc_tsv: Path,
    qc_json: Path | None = None,
    manual_candidates_tsv: Path | None = None,
    organelle: str,
    tool: str,
    phase: str,
    pga_warnings=None,
    default_transl_table="1",
):
    target_features = parse_cds_features(annotation, default_transl_table)
    references_by_gene = reference_cds_by_gene(reference_dir, default_transl_table)
    qc_rows, reference_choice_by_feature = make_qc_rows(
        target_features,
        references_by_gene,
        pga_warnings,
        organelle=organelle,
        tool=tool,
        phase=phase,
    )
    write_qc_tsv(qc_rows, qc_tsv)
    if manual_candidates_tsv is not None:
        write_manual_candidates_tsv(qc_rows, manual_candidates_tsv)
    result = {
        "annotation": str(annotation),
        "reference_dir": str(reference_dir),
        "organelle": organelle,
        "tool": tool,
        "phase": phase,
        "qc_tsv": str(qc_tsv),
        "manual_candidates_tsv": (
            str(manual_candidates_tsv) if manual_candidates_tsv else None
        ),
        "target_cds_count": len(target_features),
        "reference_gene_count": len(references_by_gene),
        "qc_row_count": len(qc_rows),
        "manual_review_candidate_count": sum(
            1 for row in qc_rows if row.manual_review_candidate
        ),
        "rows": [asdict(row) for row in qc_rows],
    }
    if qc_json is not None:
        result["qc_json"] = str(qc_json)
        write_qc_json(qc_json, result)
    return {
        **result,
        "target_features": target_features,
        "qc_rows": qc_rows,
        "reference_choice_by_feature": reference_choice_by_feature,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare organelle CDS annotations with reference CDS/proteins."
    )
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--qc-tsv", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path)
    parser.add_argument("--manual-candidates-tsv", type=Path)
    parser.add_argument("--organelle", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--phase", required=True)
    parser.add_argument("--default-transl-table", default="1")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.reference_dir.is_dir():
        raise RuntimeError(f"Reference directory does not exist: {args.reference_dir}")
    run_reference_cds_qc(
        annotation=args.annotation,
        reference_dir=args.reference_dir,
        qc_tsv=args.qc_tsv,
        qc_json=args.qc_json,
        manual_candidates_tsv=args.manual_candidates_tsv,
        organelle=args.organelle,
        tool=args.tool,
        phase=args.phase,
        default_transl_table=args.default_transl_table,
    )


if __name__ == "__main__":
    main()
