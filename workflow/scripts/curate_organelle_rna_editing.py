#!/usr/bin/env python3
import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from types import SimpleNamespace

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from organelle_annotation_utils import (  # noqa: E402
    FEATURE_INDENT,
    GENETIC_CODES,
    GENETIC_CODE_START_CODONS,
    GENETIC_CODE_STOP_CODONS,
    QUALIFIER_INDENT,
    clean_genbank_location_position,
    extract_genbank_location_sequence,
    first_genbank_qualifier_value,
    format_genbank_feature_location_lines,
    genbank_feature_has_qualifier,
    genbank_qualifier_block_value,
    genbank_qualifier_name,
    genbank_qualifier_values,
    genbank_record_locus_metadata,
    genbank_record_origin_sequence,
    parse_genbank_feature_blocks,
    normalize_genbank_cds_qualifier_order,
    split_genbank_location_arguments,
    split_genbank_records,
    translate_complete_codons,
    unwrap_genbank_location_function,
    validate_complete_cds_translation,
)
from reference_cds_qc import (  # noqa: E402
    best_reference_feature,
    protein_match_metrics,
    reference_cds_by_gene,
    reference_cds_by_product,
)


BASES = ("A", "C", "G", "T")
RNA_EDITING_EXCEPTION = "RNA editing"
RNA_EDITING_INFERENCE = "similar to RNA sequence, mRNA (same species)"
RNA_EDITING_SITE_NOTE = "C to U RNA editing"
OBSOLETE_REFERENCE_INFERRED_RNA_EDITING_SITE_NOTES = {
    "C to U RNA editing inferred from closely related species"
}
REFERENCE_SITE_EVIDENCE_NOTE_PREFIXES = (
    "inferred from reference cds evidence:",
    "inferred from reference protein evidence:",
    "inferred from reference cds:",
    "inferred from reference protein:",
)
RNA_EDITING_START_GAIN_NOTE = "start codon is created by C to U RNA editing"
RNA_EDITING_STOP_GAIN_NOTE = "stop codon is created by C to U RNA editing"
START_CODON_UNDETERMINED_NOTE = "start codon is not determined"
STOP_CODON_UNDETERMINED_NOTE = "stop codon is not determined"
TERMINAL_UNDETERMINED_NOTES = {
    START_CODON_UNDETERMINED_NOTE,
    STOP_CODON_UNDETERMINED_NOTE,
}
RNA_EDITING_OBSOLETE_NOTE_REPLACEMENTS = {
    RNA_EDITING_START_GAIN_NOTE: {
        START_CODON_UNDETERMINED_NOTE,
    },
    RNA_EDITING_STOP_GAIN_NOTE: {
        STOP_CODON_UNDETERMINED_NOTE,
    },
}
REFSEQ_ACCESSION_PREFIXES = (
    "AC_",
    "NC_",
    "NG_",
    "NM_",
    "NP_",
    "NR_",
    "NT_",
    "NW_",
    "NZ_",
    "WP_",
    "XM_",
    "XP_",
    "XR_",
    "YP_",
    "ZP_",
)
REFERENCE_INFERRED_MIN_PROTEIN_IDENTITY = 0.70
REFERENCE_INFERRED_MIN_PROTEIN_COVERAGE = 0.70
REFERENCE_PRODUCT_FILL_MIN_PROTEIN_IDENTITY = 0.70
REFERENCE_PRODUCT_FILL_MIN_PROTEIN_COVERAGE = 0.70
EVIDENCE_FIELDS = [
    "record",
    "organelle",
    "tool",
    "gene",
    "product",
    "feature_index",
    "location",
    "strand",
    "genomic_pos",
    "cds_pos",
    "coding_pos",
    "codon_index",
    "codon_pos",
    "genomic_base",
    "transcript_base",
    "edited_base",
    "genomic_codon",
    "edited_codon",
    "genomic_aa",
    "edited_aa",
    "effect",
    "candidate_source",
    "rna_depth",
    "rna_edited_reads",
    "rna_edit_fraction",
    "rna_mean_baseq",
    "rna_mean_mapq",
    "rna_sample_count_supporting",
    "dna_depth",
    "dna_alt_reads",
    "dna_alt_fraction",
    "decision",
    "original_decision",
    "rescue_reason",
    "reference_id",
    "reference_gene",
    "reference_product",
    "reference_protein_id",
    "reference_inference",
    "protein_identity",
    "protein_coverage",
    "applied",
    "note",
]
SUMMARY_FIELDS = [
    "record",
    "organelle",
    "tool",
    "gene",
    "product",
    "product_filled",
    "product_fill_reference_id",
    "product_fill_reference_gene",
    "product_fill_reference_product",
    "product_fill_reference_protein_id",
    "product_fill_protein_identity",
    "product_fill_protein_coverage",
    "feature_index",
    "location",
    "candidate_c_sites",
    "accepted_sites",
    "applied_sites",
    "rescued_sites",
    "reference_inferred_sites",
    "reference_inferred_rescue_attempted",
    "reference_inferred_candidate_sites",
    "likely_genomic_variant_sites",
    "mean_candidate_rna_depth",
    "candidate_sites_with_min_depth",
    "translation_added",
    "exception_added",
    "translation_status",
    "translation_errors",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Use RNA-seq and genomic read evidence to curate organelle CDS "
            "RNA editing exceptions and /translation qualifiers."
        )
    )
    parser.add_argument("--input-gbk", type=Path, required=True)
    parser.add_argument("--output-gbk", type=Path, required=True)
    parser.add_argument("--evidence-tsv", type=Path, required=True)
    parser.add_argument("--decisions-json", type=Path, required=True)
    parser.add_argument("--summary-tsv", type=Path, required=True)
    parser.add_argument("--input-post-curation", type=Path)
    parser.add_argument("--output-post-curation", type=Path)
    parser.add_argument("--organelle", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--transl-table", default="1")
    parser.add_argument("--rna-bam", nargs="+", type=Path, required=True)
    parser.add_argument("--dna-bam", type=Path)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--min-rna-depth", type=int, default=10)
    parser.add_argument("--min-edited-reads", type=int, default=3)
    parser.add_argument("--min-edit-fraction", type=float, default=0.10)
    parser.add_argument("--essential-rescue-min-rna-depth", type=int, default=1)
    parser.add_argument("--essential-rescue-min-edited-reads", type=int, default=1)
    parser.add_argument("--essential-rescue-min-edit-fraction", type=float, default=0.0)
    parser.add_argument(
        "--essential-rescue-max-dna-alt-fraction",
        type=float,
        default=0.10,
    )
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=30)
    parser.add_argument("--min-dna-depth", type=int, default=10)
    parser.add_argument("--max-dna-alt-fraction", type=float, default=0.10)
    return parser.parse_args()


def extract_location_positions(location: str):
    normalized = "".join(location.split())
    function_name, inner = unwrap_genbank_location_function(normalized)
    if function_name == "complement":
        return list(reversed(extract_location_positions(inner)))
    if function_name in {"join", "order"}:
        positions = []
        for part in split_genbank_location_arguments(inner):
            positions.extend(extract_location_positions(part))
        return positions
    if "," in normalized:
        positions = []
        for part in split_genbank_location_arguments(normalized):
            positions.extend(extract_location_positions(part))
        return positions

    if ":" in normalized or "^" in normalized:
        raise ValueError(f"Unsupported GenBank location: {location}")

    import re

    range_match = re.fullmatch(r"([<>]?\d+)\.\.([<>]?\d+)", normalized)
    if range_match:
        start = clean_genbank_location_position(range_match.group(1))
        end = clean_genbank_location_position(range_match.group(2))
        if start < 1 or end < 1 or start > end:
            raise ValueError(f"Unsupported GenBank location range: {location}")
        return list(range(start, end + 1))

    position_match = re.fullmatch(r"[<>]?(\d+)", normalized)
    if position_match:
        position = int(position_match.group(1))
        if position < 1:
            raise ValueError(f"Unsupported GenBank location position: {location}")
        return [position]

    raise ValueError(f"Unsupported GenBank location: {location}")


def infer_location_strand(location: str):
    normalized = "".join(location.split())
    if normalized.startswith("complement("):
        return "-"
    if "complement(" in normalized:
        return "mixed"
    return "+"


def open_bams(paths):
    import pysam

    return [pysam.AlignmentFile(str(path), "rb") for path in paths]


def close_bams(bams):
    for bam in bams:
        bam.close()


def pileup_base_counts(bams, contig, pos0, *, min_base_quality, min_mapping_quality):
    counts = {base: 0 for base in BASES}
    base_qualities = []
    mapping_qualities = []
    sample_depths = []
    sample_counts = []
    for bam in bams:
        if contig not in bam.references:
            sample_depths.append(0)
            sample_counts.append({base: 0 for base in BASES})
            continue
        sample_depth = 0
        per_sample_counts = {base: 0 for base in BASES}
        for column in bam.pileup(
            contig,
            pos0,
            pos0 + 1,
            truncate=True,
            stepper="all",
            max_depth=1_000_000,
        ):
            if column.reference_pos != pos0:
                continue
            for pileup_read in column.pileups:
                if (
                    pileup_read.is_del
                    or pileup_read.is_refskip
                    or pileup_read.query_position is None
                ):
                    continue
                alignment = pileup_read.alignment
                if (
                    alignment.is_unmapped
                    or alignment.is_secondary
                    or alignment.is_supplementary
                    or alignment.is_qcfail
                ):
                    continue
                if alignment.mapping_quality < min_mapping_quality:
                    continue
                query_position = pileup_read.query_position
                base = alignment.query_sequence[query_position].upper()
                if base not in counts:
                    continue
                if alignment.query_qualities is None:
                    base_quality = 0
                else:
                    base_quality = alignment.query_qualities[query_position]
                if base_quality < min_base_quality:
                    continue
                counts[base] += 1
                per_sample_counts[base] += 1
                sample_depth += 1
                base_qualities.append(base_quality)
                mapping_qualities.append(alignment.mapping_quality)
        sample_depths.append(sample_depth)
        sample_counts.append(per_sample_counts)
    depth = sum(counts.values())
    return {
        "counts": counts,
        "depth": depth,
        "mean_baseq": statistics.fmean(base_qualities) if base_qualities else 0.0,
        "mean_mapq": statistics.fmean(mapping_qualities) if mapping_qualities else 0.0,
        "sample_depths": sample_depths,
        "sample_counts": sample_counts,
    }


class PileupBatchCounter:
    def __init__(self, bams, *, min_base_quality, min_mapping_quality):
        self.bams = bams
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.reference_sets = [set(bam.references) for bam in bams]
        self.cache = {}

    def _new_position_state(self):
        return {
            "counts": empty_counts(),
            "baseq_sum": 0,
            "mapq_sum": 0,
            "quality_count": 0,
            "sample_depths": [0 for _bam in self.bams],
            "sample_counts": [empty_counts() for _bam in self.bams],
        }

    def _windows(self, positions0, *, full_span):
        positions0 = sorted(positions0)
        if not positions0:
            return []
        if full_span:
            return [(positions0[0], positions0[-1])]

        windows = []
        start = positions0[0]
        end = start
        max_gap = 100
        for position in positions0[1:]:
            if position - end <= max_gap:
                end = position
                continue
            windows.append((start, end))
            start = position
            end = position
        windows.append((start, end))
        return windows

    def preload(self, contig, positions0, *, full_span=False):
        positions0 = sorted({int(position) for position in positions0})
        if not positions0:
            return
        contig_cache = self.cache.setdefault(contig, {})
        missing = [position for position in positions0 if position not in contig_cache]
        if not missing:
            return
        for position in missing:
            contig_cache[position] = self._new_position_state()

        requested = set(missing)
        for bam_index, bam in enumerate(self.bams):
            if contig not in self.reference_sets[bam_index]:
                continue
            for start, end in self._windows(missing, full_span=full_span):
                for column in bam.pileup(
                    contig,
                    start,
                    end + 1,
                    truncate=True,
                    stepper="all",
                    max_depth=1_000_000,
                ):
                    if column.reference_pos not in requested:
                        continue
                    state = contig_cache[column.reference_pos]
                    per_sample_counts = state["sample_counts"][bam_index]
                    for pileup_read in column.pileups:
                        if (
                            pileup_read.is_del
                            or pileup_read.is_refskip
                            or pileup_read.query_position is None
                        ):
                            continue
                        alignment = pileup_read.alignment
                        if (
                            alignment.is_unmapped
                            or alignment.is_secondary
                            or alignment.is_supplementary
                            or alignment.is_qcfail
                        ):
                            continue
                        if alignment.mapping_quality < self.min_mapping_quality:
                            continue
                        query_position = pileup_read.query_position
                        base = alignment.query_sequence[query_position].upper()
                        if base not in state["counts"]:
                            continue
                        if alignment.query_qualities is None:
                            base_quality = 0
                        else:
                            base_quality = alignment.query_qualities[query_position]
                        if base_quality < self.min_base_quality:
                            continue
                        state["counts"][base] += 1
                        per_sample_counts[base] += 1
                        state["sample_depths"][bam_index] += 1
                        state["baseq_sum"] += base_quality
                        state["mapq_sum"] += alignment.mapping_quality
                        state["quality_count"] += 1

    def counts(self, contig, pos0):
        if pos0 not in self.cache.get(contig, {}):
            self.preload(contig, [pos0])
        state = self.cache[contig][pos0]
        depth = sum(state["counts"].values())
        quality_count = state["quality_count"]
        return {
            "counts": dict(state["counts"]),
            "depth": depth,
            "mean_baseq": (
                state["baseq_sum"] / quality_count if quality_count else 0.0
            ),
            "mean_mapq": (
                state["mapq_sum"] / quality_count if quality_count else 0.0
            ),
            "sample_depths": list(state["sample_depths"]),
            "sample_counts": [dict(counts) for counts in state["sample_counts"]],
        }


def rna_meets_accept_candidate_thresholds(
    *,
    rna_depth,
    edited_reads,
    edit_fraction,
    rna_mean_baseq,
    rna_mean_mapq,
    args,
):
    return (
        rna_depth >= args.min_rna_depth
        and edited_reads >= args.min_edited_reads
        and edit_fraction >= args.min_edit_fraction
        and rna_mean_baseq >= args.min_base_quality
        and rna_mean_mapq >= args.min_mapping_quality
    )


def classify_rna_only_decision(
    *,
    rna_depth,
    edited_reads,
    edit_fraction,
    rna_mean_baseq,
    rna_mean_mapq,
    args,
):
    if rna_depth == 0:
        return "insufficient_coverage"
    if rna_meets_accept_candidate_thresholds(
        rna_depth=rna_depth,
        edited_reads=edited_reads,
        edit_fraction=edit_fraction,
        rna_mean_baseq=rna_mean_baseq,
        rna_mean_mapq=rna_mean_mapq,
        args=args,
    ):
        return "pending_dna"
    return "reject"


def classify_decision(
    *,
    rna_depth,
    edited_reads,
    edit_fraction,
    rna_mean_baseq,
    rna_mean_mapq,
    dna_depth,
    dna_alt_fraction,
    args,
):
    dna_supported = (
        dna_depth >= args.min_dna_depth
        and dna_alt_fraction <= args.max_dna_alt_fraction
    )
    accepted = (
        rna_depth >= args.min_rna_depth
        and edited_reads >= args.min_edited_reads
        and edit_fraction >= args.min_edit_fraction
        and rna_mean_baseq >= args.min_base_quality
        and rna_mean_mapq >= args.min_mapping_quality
        and dna_supported
    )
    if accepted:
        return "accept"

    if (
        dna_depth >= args.min_dna_depth
        and dna_alt_fraction > args.max_dna_alt_fraction
    ):
        return "likely_genomic_variant"
    if rna_depth == 0:
        return "insufficient_coverage"
    if dna_depth < args.min_dna_depth:
        return "insufficient_dna_coverage"
    return "reject"


def has_essential_rescue_support(site, args):
    return (
        site["rna_depth"] >= args.essential_rescue_min_rna_depth
        and site["rna_edited_reads"] >= args.essential_rescue_min_edited_reads
        and site["rna_edit_fraction"] >= args.essential_rescue_min_edit_fraction
        and site["dna_depth"] is not None
        and site["dna_depth"] >= args.min_dna_depth
        and site["dna_alt_fraction"] <= args.essential_rescue_max_dna_alt_fraction
    )


def internal_stop_positions_from_error(error):
    prefix = "internal_stop_codon:"
    if not error.startswith(prefix):
        return set()
    positions = set()
    for text in error.removeprefix(prefix).split(","):
        try:
            positions.add(int(text))
        except ValueError:
            continue
    return positions


def essential_rescue_sites(validation_errors, site_records, terminal_codon_index, args):
    rescue_sites = []
    used_indices = set()
    internal_stop_positions = set()
    needs_start_rescue = False
    needs_terminal_stop_rescue = False

    for error in validation_errors:
        if error.startswith("invalid_start_codon:"):
            needs_start_rescue = True
        elif error.startswith("missing_terminal_stop:"):
            needs_terminal_stop_rescue = True
        elif error.startswith("internal_stop_codon:"):
            internal_stop_positions.update(internal_stop_positions_from_error(error))

    for site in site_records:
        rescue_reason = None
        if (
            needs_start_rescue
            and site["codon_index"] == 1
            and site["effect"] == "start_gain"
        ):
            rescue_reason = "invalid_start_codon"
        elif (
            needs_terminal_stop_rescue
            and site["codon_index"] == terminal_codon_index
            and site["effect"] == "stop_gain"
        ):
            rescue_reason = "missing_terminal_stop"
        elif (
            site["codon_index"] in internal_stop_positions
            and site["effect"] == "premature_stop_rescue"
        ):
            rescue_reason = "internal_stop_codon"

        if rescue_reason is None or site["cds_index"] in used_indices:
            continue
        if not has_essential_rescue_support(site, args):
            continue
        rescue_sites.append((site, rescue_reason))
        used_indices.add(site["cds_index"])

    return rescue_sites


def site_was_applied(site):
    return bool(site.get("evidence_row", {}).get("applied"))


def has_reference_inferred_support(site):
    return not site_was_applied(site)


def single_reference_rescue_site(choices):
    unique = {}
    for site in choices:
        unique[site["cds_index"]] = site
    if len(unique) != 1:
        return None
    return next(iter(unique.values()))


def reference_inferred_candidate_sites(
    validation_errors,
    site_records,
    terminal_codon_index,
):
    candidates = []
    used_indices = set()
    internal_stop_positions = set()

    for error in validation_errors:
        if error.startswith("invalid_start_codon:"):
            choices = [
                site
                for site in site_records
                if site["codon_index"] == 1
                and site["effect"] == "start_gain"
                and has_reference_inferred_support(site)
            ]
            site = single_reference_rescue_site(choices)
            if site is None:
                return []
            candidates.append((site, "invalid_start_codon"))
            used_indices.add(site["cds_index"])
        elif error.startswith("missing_terminal_stop:"):
            choices = [
                site
                for site in site_records
                if site["codon_index"] == terminal_codon_index
                and site["effect"] == "stop_gain"
                and has_reference_inferred_support(site)
            ]
            site = single_reference_rescue_site(choices)
            if site is None:
                return []
            if site["cds_index"] not in used_indices:
                candidates.append((site, "missing_terminal_stop"))
                used_indices.add(site["cds_index"])
        elif error.startswith("internal_stop_codon:"):
            internal_stop_positions.update(internal_stop_positions_from_error(error))
        elif error not in {"length_not_multiple_of_three"}:
            return []

    for codon_index in sorted(internal_stop_positions):
        choices = [
            site
            for site in site_records
            if site["codon_index"] == codon_index
            and site["effect"] == "premature_stop_rescue"
            and has_reference_inferred_support(site)
        ]
        site = single_reference_rescue_site(choices)
        if site is None:
            return []
        if site["cds_index"] not in used_indices:
            candidates.append((site, "internal_stop_codon"))
            used_indices.add(site["cds_index"])

    return candidates


def reference_inference_value(reference):
    protein_id = reference.get("protein_id", "")
    if protein_id:
        database = reference_sequence_database(protein_id)
        return f"similar to AA sequence:{database}:{protein_id}"
    accession = reference.get("accession", "") or reference.get("record_id", "")
    if accession:
        database = reference_sequence_database(accession)
        return f"similar to AA sequence:{database}:{accession}"
    return "similar to AA sequence"


def reference_sequence_database(accession):
    if accession.startswith(REFSEQ_ACCESSION_PREFIXES):
        return "RefSeq"
    return "INSD"


def reference_evidence_accession(reference):
    return (
        reference.get("protein_id", "")
        or reference.get("accession", "")
        or reference.get("record_id", "")
    )


def reference_evidence_basis(reference):
    accession = reference_evidence_accession(reference)
    if not accession:
        return ""
    database = reference_sequence_database(accession)
    return f"{database}:{accession}"


def reference_site_evidence_note(reference):
    reference = reference or {}
    evidence_basis = reference_evidence_basis(reference)
    evidence_kind = "protein" if reference.get("protein_id", "") else "CDS"
    if evidence_basis:
        return f"inferred from reference {evidence_kind}: {evidence_basis}"
    return "inferred from reference CDS"


def best_reference_for_feature(
    *,
    gene,
    product,
    coding_sequence,
    cds_length,
    table_id,
    references_by_gene,
    references_by_product,
):
    if not references_by_gene:
        return None
    query_protein = translate_cds_for_genbank(coding_sequence, table_id)
    target = SimpleNamespace(
        gene=gene,
        product=product,
        protein=query_protein,
        length=cds_length,
    )
    return best_reference_feature(
        target,
        references_by_gene,
        references_by_product,
    )


def product_qualifier_needs_fill(product_values):
    return not product_values or all(not str(value).strip() for value in product_values)


def reference_product_fill_candidate(
    *,
    gene,
    product_values,
    coding_sequence,
    cds_length,
    table_id,
    references_by_gene,
):
    if not references_by_gene or not product_qualifier_needs_fill(product_values):
        return None
    if not gene or gene == "?":
        return None

    query_protein = translate_cds_for_genbank(coding_sequence, table_id)
    if not query_protein:
        return None

    candidates = []
    for reference in references_by_gene.get(gene, []):
        reference_product = str(reference.get("product", "")).strip()
        reference_protein = reference.get("protein", "")
        if not reference_product or not reference_protein:
            continue
        metrics = protein_match_metrics(query_protein, reference_protein)
        if (
            metrics["protein_identity"] is None
            or metrics["protein_coverage"] is None
            or metrics["protein_identity"] < REFERENCE_PRODUCT_FILL_MIN_PROTEIN_IDENTITY
            or metrics["protein_coverage"] < REFERENCE_PRODUCT_FILL_MIN_PROTEIN_COVERAGE
        ):
            continue
        reference_length = reference.get("length") or 0
        length_delta = abs(reference_length - cds_length) if reference_length else 0
        candidates.append(
            {
                "reference": reference,
                "product": reference_product,
                "metrics": metrics,
                "score": (
                    metrics["protein_identity"],
                    metrics["protein_coverage"],
                    -length_delta,
                ),
            }
        )

    if not candidates:
        return None
    candidates.sort(key=lambda candidate: candidate["score"], reverse=True)
    best_score = candidates[0]["score"]
    best_products = {
        candidate["product"].lower()
        for candidate in candidates
        if candidate["score"] == best_score
    }
    if len(best_products) > 1:
        return None
    return candidates[0]


def reference_product_fill_from_reference(
    *,
    gene,
    product_values,
    reference,
    metrics,
):
    if not product_qualifier_needs_fill(product_values):
        return None
    if reference is None or reference.get("gene") != gene:
        return None
    reference_product = str(reference.get("product", "")).strip()
    if not reference_product:
        return None
    if (
        metrics.get("protein_identity") is None
        or metrics.get("protein_coverage") is None
        or metrics["protein_identity"] < REFERENCE_PRODUCT_FILL_MIN_PROTEIN_IDENTITY
        or metrics["protein_coverage"] < REFERENCE_PRODUCT_FILL_MIN_PROTEIN_COVERAGE
    ):
        return None
    return {
        "reference": reference,
        "product": reference_product,
        "metrics": metrics,
    }


def reference_inferred_rescue(
    *,
    validation_errors,
    site_records,
    edited_cds,
    codon_offset,
    terminal_codon_index,
    table_id,
    gene,
    product,
    references_by_gene,
    references_by_product,
    rescue_sites=None,
):
    if rescue_sites is None:
        rescue_sites = reference_inferred_candidate_sites(
            validation_errors,
            site_records,
            terminal_codon_index,
        )
    if not rescue_sites:
        return None

    rescued_cds = list(edited_cds)
    for site, _reason in rescue_sites:
        rescued_cds[site["cds_index"]] = "T"
    rescued_coding_sequence = "".join(rescued_cds)[codon_offset:]
    rescued_validation = validate_complete_cds_translation(
        rescued_coding_sequence,
        table_id,
    )
    if not rescued_validation["valid"]:
        return None

    reference = best_reference_for_feature(
        gene=gene,
        product=product,
        coding_sequence=rescued_coding_sequence,
        cds_length=len(rescued_cds),
        table_id=table_id,
        references_by_gene=references_by_gene,
        references_by_product=references_by_product,
    )
    if reference is None or not reference.get("protein"):
        return None

    translation = translate_cds_for_genbank(rescued_coding_sequence, table_id)
    metrics = protein_match_metrics(translation, reference.get("protein", ""))
    if (
        metrics["protein_identity"] is None
        or metrics["protein_coverage"] is None
        or metrics["protein_identity"] < REFERENCE_INFERRED_MIN_PROTEIN_IDENTITY
        or metrics["protein_coverage"] < REFERENCE_INFERRED_MIN_PROTEIN_COVERAGE
    ):
        return None

    return {
        "sites": rescue_sites,
        "coding_sequence": rescued_coding_sequence,
        "translation": translation,
        "reference": reference,
        "metrics": metrics,
        "inference": reference_inference_value(reference),
    }


def codon_effect(
    *,
    codon_index,
    genomic_codon,
    edited_codon,
    table_id,
):
    genomic_aa = GENETIC_CODES[table_id].get(genomic_codon, "X")
    edited_aa = GENETIC_CODES[table_id].get(edited_codon, "X")
    stop_codons = GENETIC_CODE_STOP_CODONS[table_id]
    start_codons = GENETIC_CODE_START_CODONS[table_id]
    if codon_index == 1 and edited_codon in start_codons and genomic_codon not in start_codons:
        effect = "start_gain"
    elif edited_codon in stop_codons and genomic_codon not in stop_codons:
        effect = "stop_gain" if genomic_codon != edited_codon else "none"
    elif genomic_codon in stop_codons and edited_codon not in stop_codons:
        effect = "premature_stop_rescue"
    elif genomic_aa == edited_aa:
        effect = "synonymous"
    else:
        effect = "nonsynonymous"
    return effect, genomic_aa, edited_aa


def translate_cds_for_genbank(sequence, table_id):
    table_id = str(table_id)
    protein, _ambiguous = translate_complete_codons(sequence, table_id)
    codons = [
        sequence[index : index + 3].upper()
        for index in range(0, len(sequence) - len(sequence) % 3, 3)
    ]
    if codons and codons[0] in GENETIC_CODE_START_CODONS[table_id]:
        protein = "M" + protein[1:]
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def validation_has_error(errors, prefix):
    return any(error == prefix or error.startswith(f"{prefix}:") for error in errors)


def validation_errors_allow_undetermined_terminal(errors):
    allowed = {
        "invalid_start_codon",
        "missing_terminal_stop",
        "length_not_multiple_of_three",
    }
    if not errors:
        return False
    if any(
        not any(error == prefix or error.startswith(f"{prefix}:") for prefix in allowed)
        for error in errors
    ):
        return False
    has_incomplete_codon = validation_has_error(errors, "length_not_multiple_of_three")
    has_missing_stop = validation_has_error(errors, "missing_terminal_stop")
    if has_incomplete_codon and not has_missing_stop:
        return False
    return validation_has_error(errors, "invalid_start_codon") or has_missing_stop


def validate_cds_translation_allowing_undetermined_terminal(sequence, table_id):
    validation = validate_complete_cds_translation(sequence, table_id)
    if validation["valid"]:
        return {
            "valid": True,
            "start_undetermined": False,
            "stop_undetermined": False,
            "errors": [],
            "warnings": validation["warnings"],
        }
    errors = validation["errors"]
    start_undetermined = validation_has_error(errors, "invalid_start_codon")
    stop_undetermined = validation_has_error(errors, "missing_terminal_stop")
    if not validation_errors_allow_undetermined_terminal(errors):
        return {
            "valid": False,
            "start_undetermined": start_undetermined,
            "stop_undetermined": stop_undetermined,
            "errors": errors,
            "warnings": validation["warnings"],
        }
    return {
        "valid": True,
        "start_undetermined": start_undetermined,
        "stop_undetermined": stop_undetermined,
        "errors": errors,
        "warnings": validation["warnings"],
    }


def add_partial_markers_to_simple_location(location, *, five_prime, three_prime):
    import re

    range_match = re.fullmatch(r"([<>]?)(\d+)\.\.([<>]?)(\d+)", location)
    if range_match:
        left_marker = "<" if five_prime else range_match.group(1)
        right_marker = ">" if three_prime else range_match.group(3)
        return (
            f"{left_marker}{range_match.group(2)}.."
            f"{right_marker}{range_match.group(4)}"
        )

    position_match = re.fullmatch(r"([<>]?)(\d+)", location)
    if position_match:
        marker = "<" if five_prime else ">" if three_prime else position_match.group(1)
        return f"{marker}{position_match.group(2)}"

    return location


def add_partial_markers_to_location(location, *, five_prime=False, three_prime=False):
    normalized = "".join(location.split())
    function_name, inner = unwrap_genbank_location_function(normalized)
    if function_name == "complement":
        return (
            "complement("
            + add_partial_markers_to_location(
                inner,
                five_prime=three_prime,
                three_prime=five_prime,
            )
            + ")"
        )
    if function_name in {"join", "order"}:
        parts = split_genbank_location_arguments(inner)
        if not parts:
            return normalized
        marked_parts = [
            add_partial_markers_to_location(
                part,
                five_prime=five_prime and index == 0,
                three_prime=three_prime and index == len(parts) - 1,
            )
            for index, part in enumerate(parts)
        ]
        return f"{function_name}(" + ",".join(marked_parts) + ")"
    if "," in normalized:
        parts = split_genbank_location_arguments(normalized)
        return ",".join(
            add_partial_markers_to_location(
                part,
                five_prime=five_prime and index == 0,
                three_prime=three_prime and index == len(parts) - 1,
            )
            for index, part in enumerate(parts)
        )
    return add_partial_markers_to_simple_location(
        normalized,
        five_prime=five_prime,
        three_prime=three_prime,
    )


def replace_cds_block_location(block_lines, location):
    qualifier_start = len(block_lines)
    for index, line in enumerate(block_lines):
        if genbank_qualifier_name(line) is not None:
            qualifier_start = index
            break
    feature_key = block_lines[0][len(FEATURE_INDENT) : len(FEATURE_INDENT) + 16].strip()
    return [
        *format_genbank_feature_location_lines(feature_key or "CDS", location),
        *block_lines[qualifier_start:],
    ]


def format_quoted_qualifier(name, value, width=58):
    prefix = f"/{name}=\""
    if len(value) <= width:
        return [f"{QUALIFIER_INDENT}{prefix}{value}\"\n"]
    lines = [f"{QUALIFIER_INDENT}{prefix}{value[:width]}\n"]
    value = value[width:]
    while len(value) > width:
        lines.append(f"{QUALIFIER_INDENT}{value[:width]}\n")
        value = value[width:]
    lines.append(f"{QUALIFIER_INDENT}{value}\"\n")
    return lines


def format_simple_qualifier(name, value):
    return [f"{QUALIFIER_INDENT}/{name}={value}\n"]


def qualifier_block_end(block_lines, start_index):
    end_index = start_index + 1
    while (
        end_index < len(block_lines)
        and genbank_qualifier_name(block_lines[end_index]) is None
    ):
        end_index += 1
    return end_index


def fill_cds_product_qualifier(block_lines, product):
    cleaned = remove_qualifiers(block_lines, {"product"})
    insert_index = len(cleaned)
    for index, line in enumerate(cleaned):
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name in {"gene", "locus_tag", "old_locus_tag"}:
            insert_index = qualifier_block_end(cleaned, index)
    return [
        *cleaned[:insert_index],
        *format_quoted_qualifier("product", product),
        *cleaned[insert_index:],
    ]


def product_fill_reference_summary(product_fill):
    reference = product_fill["reference"]
    metrics = product_fill["metrics"]
    return {
        "product_filled": True,
        "product_fill_reference_id": reference.get("record_id", ""),
        "product_fill_reference_gene": reference.get("gene", ""),
        "product_fill_reference_product": reference.get("product", ""),
        "product_fill_reference_protein_id": reference.get("protein_id", ""),
        "product_fill_protein_identity": f"{metrics['protein_identity']:.6f}",
        "product_fill_protein_coverage": f"{metrics['protein_coverage']:.6f}",
    }


def product_fill_reference_decision(product_fill):
    reference = product_fill["reference"]
    metrics = product_fill["metrics"]
    return {
        "reference_id": reference.get("record_id", ""),
        "reference_accession": reference.get("accession", ""),
        "reference_gene": reference.get("gene", ""),
        "reference_product": reference.get("product", ""),
        "reference_protein_id": reference.get("protein_id", ""),
        "reference_path": reference.get("path", ""),
        "protein_identity": metrics["protein_identity"],
        "protein_coverage": metrics["protein_coverage"],
    }


def remove_qualifiers(block_lines, qualifiers):
    qualifier_set = set(qualifiers)
    cleaned = []
    skipping = False
    for line in block_lines:
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is not None:
            if qualifier_name in qualifier_set:
                skipping = True
                continue
            skipping = False
            cleaned.append(line)
            continue
        if skipping:
            continue
        cleaned.append(line)
    return cleaned


def qualifier_value_from_line(line, qualifier):
    return genbank_qualifier_block_value(qualifier, [line])


def remove_qualifier_values(block_lines, qualifier, values):
    value_set = {value.lower() for value in values}
    cleaned = []
    index = 0
    while index < len(block_lines):
        line = block_lines[index]
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is None:
            cleaned.append(line)
            index += 1
            continue

        next_index = index + 1
        while (
            next_index < len(block_lines)
            and genbank_qualifier_name(block_lines[next_index]) is None
        ):
            next_index += 1
        qualifier_lines = block_lines[index:next_index]
        if qualifier_name == qualifier:
            value = genbank_qualifier_block_value(qualifier, qualifier_lines)
            if value is not None and value.lower() in value_set:
                index = next_index
                continue
        cleaned.extend(qualifier_lines)
        index = next_index
    return cleaned


def remove_stale_reference_site_notes(block_lines, desired_notes):
    desired_note_set = {note.lower() for note in desired_notes}
    cleaned = []
    index = 0
    while index < len(block_lines):
        line = block_lines[index]
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is None:
            cleaned.append(line)
            index += 1
            continue

        next_index = index + 1
        while (
            next_index < len(block_lines)
            and genbank_qualifier_name(block_lines[next_index]) is None
        ):
            next_index += 1
        qualifier_lines = block_lines[index:next_index]
        if qualifier_name == "note":
            value = genbank_qualifier_block_value("note", qualifier_lines)
            if value is not None:
                value_lower = value.lower()
                if (
                    any(
                        value_lower.startswith(prefix)
                        for prefix in REFERENCE_SITE_EVIDENCE_NOTE_PREFIXES
                    )
                    and value_lower not in desired_note_set
                ):
                    index = next_index
                    continue
        cleaned.extend(qualifier_lines)
        index = next_index
    return cleaned


def listify_qualifier_values(values):
    if values is None:
        return []
    if isinstance(values, str):
        return [values]
    return list(values)


def update_cds_block_lines(
    block_lines,
    *,
    translation,
    table_id,
    add_rna_editing_exception,
    location=None,
    recalculate_terminal_notes=False,
    rna_editing_note=None,
    rna_editing_inferences=None,
    extra_notes=None,
):
    cleaned = remove_qualifiers(block_lines, {"translation"})
    if location is not None:
        cleaned = replace_cds_block_location(cleaned, location)
    if recalculate_terminal_notes:
        cleaned = remove_qualifier_values(
            cleaned,
            "note",
            TERMINAL_UNDETERMINED_NOTES,
        )
    existing_exceptions = genbank_qualifier_values(cleaned, "exception")
    has_rna_editing_exception = any(
        value.lower() == RNA_EDITING_EXCEPTION.lower() for value in existing_exceptions
    )
    if add_rna_editing_exception and not has_rna_editing_exception:
        cleaned.extend(format_quoted_qualifier("exception", RNA_EDITING_EXCEPTION))
    if add_rna_editing_exception:
        rna_editing_notes = listify_qualifier_values(rna_editing_note)
        obsolete_notes = set()
        for note in rna_editing_notes:
            obsolete_notes.update(RNA_EDITING_OBSOLETE_NOTE_REPLACEMENTS.get(note, ()))
        if obsolete_notes:
            cleaned = remove_qualifier_values(cleaned, "note", obsolete_notes)
        existing_notes = {
            value.lower() for value in genbank_qualifier_values(cleaned, "note")
        }
        for note in rna_editing_notes:
            if note.lower() not in existing_notes:
                cleaned.extend(format_quoted_qualifier("note", note))
                existing_notes.add(note.lower())
        existing_inferences = {
            value.lower()
            for value in genbank_qualifier_values(cleaned, "inference")
        }
        for inference in rna_editing_inferences or [RNA_EDITING_INFERENCE]:
            if inference.lower() not in existing_inferences:
                cleaned.extend(format_quoted_qualifier("inference", inference))
                existing_inferences.add(inference.lower())
    extra_notes = listify_qualifier_values(extra_notes)
    if extra_notes:
        existing_notes = {
            value.lower() for value in genbank_qualifier_values(cleaned, "note")
        }
        for note in extra_notes:
            if note.lower() not in existing_notes:
                cleaned.extend(format_quoted_qualifier("note", note))
                existing_notes.add(note.lower())
    if table_id != "1" and not any(
        genbank_qualifier_name(line) == "transl_table" for line in cleaned
    ):
        cleaned.extend(format_simple_qualifier("transl_table", table_id))
    cleaned.extend(format_quoted_qualifier("translation", translation))
    return cleaned


def rna_editing_cds_note(editing_effects, applied_decisions):
    _ = applied_decisions
    notes = []
    if "start_gain" in editing_effects:
        notes.append(RNA_EDITING_START_GAIN_NOTE)
    if "stop_gain" in editing_effects:
        notes.append(RNA_EDITING_STOP_GAIN_NOTE)
    if not notes:
        return None
    if len(notes) == 1:
        return notes[0]
    return notes


def undetermined_terminal_translation_status(start_undetermined, stop_undetermined):
    if start_undetermined and stop_undetermined:
        return "updated_with_undetermined_start_and_stop"
    if start_undetermined:
        return "updated_with_undetermined_start"
    return "updated_with_undetermined_stop"


def should_recalculate_terminal_notes(args):
    return args.organelle == "mitochondrion" and args.tool == "pmga"


def rna_editing_site_location(site):
    position = site["genomic_pos"]
    if site.get("genomic_base") == "G":
        return f"complement({position})"
    return str(position)


def rna_editing_site_notes(site):
    notes = site.get("site_notes", site.get("site_note", RNA_EDITING_SITE_NOTE))
    notes = listify_qualifier_values(notes)
    return notes or [RNA_EDITING_SITE_NOTE]


def rna_editing_site_feature_lines(site):
    location = rna_editing_site_location(site)
    return [
        f"{FEATURE_INDENT}{'misc_feature':<16}{location}\n",
        *[
            line
            for note in rna_editing_site_notes(site)
            for line in format_quoted_qualifier("note", note)
        ],
    ]


def existing_rna_editing_site_feature_registry(blocks):
    existing = set()
    block_lines_by_key = {}
    accepted_notes = {RNA_EDITING_SITE_NOTE.lower()}
    accepted_notes.update(
        note.lower() for note in OBSOLETE_REFERENCE_INFERRED_RNA_EDITING_SITE_NOTES
    )
    for block in blocks:
        if block.key != "misc_feature":
            continue
        notes = {
            value.lower() for value in genbank_qualifier_values(block.lines, "note")
        }
        for note in notes & accepted_notes:
            key = (block.location, RNA_EDITING_SITE_NOTE.lower())
            existing.add(key)
            block_lines_by_key.setdefault(key, block.lines)
    return existing, block_lines_by_key


def existing_rna_editing_site_features(blocks):
    existing, _block_lines_by_key = existing_rna_editing_site_feature_registry(blocks)
    return existing


def merge_rna_editing_site_feature_notes(block_lines, desired_notes):
    cleaned = remove_qualifier_values(
        block_lines,
        "note",
        OBSOLETE_REFERENCE_INFERRED_RNA_EDITING_SITE_NOTES,
    )
    cleaned = remove_stale_reference_site_notes(cleaned, desired_notes)
    if cleaned != block_lines:
        block_lines[:] = cleaned
    existing_notes = {
        value.lower() for value in genbank_qualifier_values(block_lines, "note")
    }
    for note in desired_notes:
        if note.lower() in existing_notes:
            continue
        block_lines.extend(format_quoted_qualifier("note", note))
        existing_notes.add(note.lower())


def new_rna_editing_site_feature_blocks(
    feature_decision,
    existing_site_features,
    existing_site_feature_lines=None,
):
    existing_site_feature_lines = existing_site_feature_lines or {}
    blocks = []
    for site in feature_decision["accepted_sites"]:
        location = rna_editing_site_location(site)
        site_notes = rna_editing_site_notes(site)
        site_note = site_notes[0]
        key = (location, site_note.lower())
        if key in existing_site_features:
            existing_lines = existing_site_feature_lines.get(key)
            if existing_lines is not None:
                merge_rna_editing_site_feature_notes(existing_lines, site_notes)
            continue
        existing_site_features.add(key)
        block_lines = rna_editing_site_feature_lines(site)
        existing_site_feature_lines[key] = block_lines
        blocks.append(block_lines)
    return blocks


def feature_candidate_sources(tool, exceptions):
    sources = ["all_CDS_scan"]
    if any(value.lower() == "rna editing" for value in exceptions):
        sources.append(f"{tool}_existing_RNA_editing_exception")
    return ";".join(sources)


def empty_counts():
    return {base: 0 for base in BASES}


def make_accepted_site(site, *, decision, rescue_reason="", reference=None, metrics=None):
    accepted_site = {
        "genomic_pos": site["genomic_pos"],
        "cds_pos": site["cds_pos"],
        "codon_index": site["codon_index"],
        "codon_pos": site["codon_pos"],
        "genomic_base": site["genomic_base"],
        "edited_base": site["edited_base"],
        "effect": site["effect"],
        "decision": decision,
        "rna_depth": site["rna_depth"],
        "rna_edited_reads": site["rna_edited_reads"],
        "rna_edit_fraction": site["rna_edit_fraction"],
        "dna_depth": site["dna_depth"],
        "dna_alt_fraction": site["dna_alt_fraction"],
    }
    if rescue_reason:
        accepted_site["rescue_reason"] = rescue_reason
    if decision == "reference_inferred":
        accepted_site["site_notes"] = [
            RNA_EDITING_SITE_NOTE,
            reference_site_evidence_note(reference),
        ]
        accepted_site["reference_id"] = reference.get("record_id", "") if reference else ""
        accepted_site["reference_gene"] = reference.get("gene", "") if reference else ""
        accepted_site["reference_product"] = (
            reference.get("product", "") if reference else ""
        )
        accepted_site["reference_protein_id"] = (
            reference.get("protein_id", "") if reference else ""
        )
        if metrics:
            accepted_site["protein_identity"] = metrics.get("protein_identity")
            accepted_site["protein_coverage"] = metrics.get("protein_coverage")
    return accepted_site


def apply_dna_counts_to_site(site, dna_counts):
    edited_aligned_base = site["edited_aligned_base"]
    dna_depth = dna_counts["depth"]
    dna_alt_reads = dna_counts["counts"][edited_aligned_base]
    dna_alt_fraction = dna_alt_reads / dna_depth if dna_depth > 0 else 1.0
    site["dna_depth"] = dna_depth
    site["dna_alt_reads"] = dna_alt_reads
    site["dna_alt_fraction"] = dna_alt_fraction
    site["dna_checked"] = True
    evidence_row = site["evidence_row"]
    evidence_row["dna_depth"] = dna_depth
    evidence_row["dna_alt_reads"] = dna_alt_reads
    evidence_row["dna_alt_fraction"] = f"{dna_alt_fraction:.6f}"


def ensure_dna_counts_for_sites(sites, *, dna_pileups, record_name):
    sites_to_load = [site for site in sites if not site.get("dna_checked")]
    if not sites_to_load:
        return
    dna_pileups.preload(
        record_name,
        [site["genomic_pos"] - 1 for site in sites_to_load],
    )
    for site in sites_to_load:
        dna_counts = dna_pileups.counts(record_name, site["genomic_pos"] - 1)
        apply_dna_counts_to_site(site, dna_counts)


def has_essential_rescue_rna_support(site, args):
    return (
        site["rna_depth"] >= args.essential_rescue_min_rna_depth
        and site["rna_edited_reads"] >= args.essential_rescue_min_edited_reads
        and site["rna_edit_fraction"] >= args.essential_rescue_min_edit_fraction
    )


def essential_rescue_dna_candidate_sites(
    validation_errors,
    site_records,
    terminal_codon_index,
    args,
):
    candidates = []
    internal_stop_positions = set()
    needs_start_rescue = False
    needs_terminal_stop_rescue = False

    for error in validation_errors:
        if error.startswith("invalid_start_codon:"):
            needs_start_rescue = True
        elif error.startswith("missing_terminal_stop:"):
            needs_terminal_stop_rescue = True
        elif error.startswith("internal_stop_codon:"):
            internal_stop_positions.update(internal_stop_positions_from_error(error))

    for site in site_records:
        if not has_essential_rescue_rna_support(site, args):
            continue
        if (
            needs_start_rescue
            and site["codon_index"] == 1
            and site["effect"] == "start_gain"
        ):
            candidates.append(site)
        elif (
            needs_terminal_stop_rescue
            and site["codon_index"] == terminal_codon_index
            and site["effect"] == "stop_gain"
        ):
            candidates.append(site)
        elif (
            site["codon_index"] in internal_stop_positions
            and site["effect"] == "premature_stop_rescue"
        ):
            candidates.append(site)

    return candidates


def candidate_positions_for_cds_block(block, record_sequence):
    is_pseudo = genbank_feature_has_qualifier(
        block.lines,
        "pseudo",
    ) or genbank_feature_has_qualifier(block.lines, "pseudogene")
    if is_pseudo:
        return []
    try:
        cds_sequence = extract_genbank_location_sequence(
            block.location,
            record_sequence,
        )
        genomic_positions = extract_location_positions(block.location)
    except ValueError:
        return []
    if len(cds_sequence) != len(genomic_positions):
        return []

    codon_start_text = first_genbank_qualifier_value(block.lines, "codon_start", "1")
    try:
        codon_offset = int(codon_start_text) - 1
    except (TypeError, ValueError):
        codon_offset = 0
    if codon_offset not in {0, 1, 2}:
        codon_offset = 0
    complete_coding_length = len(cds_sequence[codon_offset:])
    complete_coding_length -= complete_coding_length % 3
    coding_end = codon_offset + complete_coding_length

    positions = []
    for cds_index in range(codon_offset, coding_end):
        if cds_sequence[cds_index].upper() != "C":
            continue
        genomic_pos = genomic_positions[cds_index]
        genomic_base = record_sequence[genomic_pos - 1].upper()
        if genomic_base in {"C", "G"}:
            positions.append(genomic_pos - 1)
    return positions


def candidate_positions_for_record(blocks, record_sequence):
    positions = []
    for block in blocks:
        if block.key == "CDS":
            positions.extend(candidate_positions_for_cds_block(block, record_sequence))
    return positions


def evaluate_cds_feature(
    *,
    record_name,
    record_sequence,
    block,
    args,
    rna_pileups,
    dna_pileups,
    references_by_gene=None,
    references_by_product=None,
):
    table_id = str(
        first_genbank_qualifier_value(block.lines, "transl_table", args.transl_table)
    )
    if table_id not in GENETIC_CODES:
        table_id = str(args.transl_table)
    gene = first_genbank_qualifier_value(block.lines, "gene", "?")
    product_values = genbank_qualifier_values(block.lines, "product")
    product = product_values[0] if product_values else ""
    exceptions = genbank_qualifier_values(block.lines, "exception")
    candidate_source = feature_candidate_sources(args.tool, exceptions)
    is_pseudo = genbank_feature_has_qualifier(
        block.lines,
        "pseudo",
    ) or genbank_feature_has_qualifier(block.lines, "pseudogene")
    strand = infer_location_strand(block.location)
    recalculate_terminal_notes = should_recalculate_terminal_notes(args)

    summary = {
        "record": record_name,
        "organelle": args.organelle,
        "tool": args.tool,
        "gene": gene,
        "product": product,
        "product_filled": False,
        "product_fill_reference_id": "",
        "product_fill_reference_gene": "",
        "product_fill_reference_product": "",
        "product_fill_reference_protein_id": "",
        "product_fill_protein_identity": "",
        "product_fill_protein_coverage": "",
        "feature_index": block.feature_index,
        "location": block.location,
        "candidate_c_sites": 0,
        "accepted_sites": 0,
        "applied_sites": 0,
        "rescued_sites": 0,
        "reference_inferred_sites": 0,
        "reference_inferred_rescue_attempted": False,
        "reference_inferred_candidate_sites": 0,
        "likely_genomic_variant_sites": 0,
        "mean_candidate_rna_depth": 0,
        "candidate_sites_with_min_depth": 0,
        "translation_added": False,
        "exception_added": False,
        "translation_status": "not_evaluated",
        "translation_errors": "",
    }
    feature_decision = {
        "record": record_name,
        "feature_index": block.feature_index,
        "gene": gene,
        "product": product,
        "product_filled": False,
        "product_fill_reference": None,
        "location": block.location,
        "candidate_source": candidate_source,
        "accepted_sites": [],
        "translation_added": False,
        "exception_added": False,
        "translation_status": "not_evaluated",
        "translation_errors": [],
    }

    evidence_rows = []
    site_records = []
    if is_pseudo:
        summary["translation_status"] = "skipped_pseudo"
        feature_decision["translation_status"] = summary["translation_status"]
        return block.lines, evidence_rows, summary, feature_decision

    try:
        cds_sequence = extract_genbank_location_sequence(
            block.location,
            record_sequence,
        )
        genomic_positions = extract_location_positions(block.location)
    except ValueError as error:
        summary["translation_status"] = "location_parse_error"
        summary["translation_errors"] = str(error)
        feature_decision["translation_status"] = summary["translation_status"]
        feature_decision["translation_errors"] = [str(error)]
        return block.lines, evidence_rows, summary, feature_decision

    if len(cds_sequence) != len(genomic_positions):
        error = (
            "CDS sequence length and parsed genomic position count differ: "
            f"{len(cds_sequence)} vs {len(genomic_positions)}"
        )
        summary["translation_status"] = "location_parse_error"
        summary["translation_errors"] = error
        feature_decision["translation_status"] = summary["translation_status"]
        feature_decision["translation_errors"] = [error]
        return block.lines, evidence_rows, summary, feature_decision

    codon_start_text = first_genbank_qualifier_value(block.lines, "codon_start", "1")
    try:
        codon_offset = int(codon_start_text) - 1
    except (TypeError, ValueError):
        codon_offset = 0
    if codon_offset not in {0, 1, 2}:
        codon_offset = 0

    block_lines = block.lines
    product_fill = reference_product_fill_candidate(
        gene=gene,
        product_values=product_values,
        coding_sequence=cds_sequence[codon_offset:],
        cds_length=len(cds_sequence),
        table_id=table_id,
        references_by_gene=references_by_gene,
    )
    if product_fill is not None:
        product = product_fill["product"]
        product_values = [product]
        block_lines = fill_cds_product_qualifier(block_lines, product)
        summary["product"] = product
        summary.update(product_fill_reference_summary(product_fill))
        feature_decision["product"] = product
        feature_decision["product_filled"] = True
        feature_decision["product_fill_reference"] = product_fill_reference_decision(
            product_fill
        )

    edited_cds = list(cds_sequence)
    rna_depths = []
    applied_indices = set()
    editing_effects = []
    applied_decisions = []
    complete_coding_length = len(cds_sequence[codon_offset:])
    complete_coding_length -= complete_coding_length % 3
    coding_end = codon_offset + complete_coding_length
    normal_dna_candidate_sites = []

    for cds_index in range(codon_offset, coding_end):
        transcript_base = cds_sequence[cds_index].upper()
        if transcript_base != "C":
            continue
        codon_index = ((cds_index - codon_offset) // 3) + 1
        codon_pos = ((cds_index - codon_offset) % 3) + 1
        codon_start = codon_offset + ((cds_index - codon_offset) // 3) * 3
        genomic_codon = cds_sequence[codon_start : codon_start + 3].upper()
        edited_codon_list = list(genomic_codon)
        edited_codon_list[codon_pos - 1] = "T"
        edited_codon = "".join(edited_codon_list)
        genomic_pos = genomic_positions[cds_index]
        genomic_base = record_sequence[genomic_pos - 1].upper()
        if genomic_base == "C":
            edited_aligned_base = "T"
        elif genomic_base == "G":
            edited_aligned_base = "A"
        else:
            continue

        effect, genomic_aa, edited_aa = codon_effect(
            codon_index=codon_index,
            genomic_codon=genomic_codon,
            edited_codon=edited_codon,
            table_id=table_id,
        )
        rna_counts = rna_pileups.counts(record_name, genomic_pos - 1)
        rna_depth = rna_counts["depth"]
        rna_edited_reads = rna_counts["counts"][edited_aligned_base]
        rna_edit_fraction = (
            rna_edited_reads / rna_depth if rna_depth > 0 else 0.0
        )
        decision = classify_rna_only_decision(
            rna_depth=rna_depth,
            edited_reads=rna_edited_reads,
            edit_fraction=rna_edit_fraction,
            rna_mean_baseq=rna_counts["mean_baseq"],
            rna_mean_mapq=rna_counts["mean_mapq"],
            args=args,
        )
        site = {
            "cds_index": cds_index,
            "genomic_pos": genomic_pos,
            "cds_pos": cds_index + 1,
            "codon_index": codon_index,
            "codon_pos": codon_pos,
            "genomic_base": genomic_base,
            "edited_base": edited_aligned_base,
            "effect": effect,
            "rna_depth": rna_depth,
            "rna_edited_reads": rna_edited_reads,
            "rna_edit_fraction": rna_edit_fraction,
            "rna_mean_baseq": rna_counts["mean_baseq"],
            "rna_mean_mapq": rna_counts["mean_mapq"],
            "dna_depth": None,
            "dna_alt_reads": None,
            "dna_alt_fraction": None,
            "dna_checked": False,
            "edited_aligned_base": edited_aligned_base,
        }

        rna_depths.append(rna_depth)
        summary["candidate_c_sites"] += 1
        summary["candidate_sites_with_min_depth"] += int(
            rna_depth >= args.min_rna_depth
        )
        note = ""
        if decision == "reject":
            note = (
                "DNA/HiFi pileup was not evaluated because RNA support did "
                "not pass normal RNA-editing candidate thresholds."
            )
        evidence_row = {
            "record": record_name,
            "organelle": args.organelle,
            "tool": args.tool,
            "gene": gene,
            "product": product,
            "feature_index": block.feature_index,
            "location": block.location,
            "strand": strand,
            "genomic_pos": genomic_pos,
            "cds_pos": cds_index + 1,
            "coding_pos": cds_index - codon_offset + 1,
            "codon_index": codon_index,
            "codon_pos": codon_pos,
            "genomic_base": genomic_base,
            "transcript_base": transcript_base,
            "edited_base": edited_aligned_base,
            "genomic_codon": genomic_codon,
            "edited_codon": edited_codon,
            "genomic_aa": genomic_aa,
            "edited_aa": edited_aa,
            "effect": effect,
            "candidate_source": candidate_source,
            "rna_depth": rna_depth,
            "rna_edited_reads": rna_edited_reads,
            "rna_edit_fraction": f"{rna_edit_fraction:.6f}",
            "rna_mean_baseq": f"{rna_counts['mean_baseq']:.3f}",
            "rna_mean_mapq": f"{rna_counts['mean_mapq']:.3f}",
            "rna_sample_count_supporting": sum(
                1
                for sample_counts in rna_counts["sample_counts"]
                if sample_counts[edited_aligned_base] > 0
            ),
            "dna_depth": "",
            "dna_alt_reads": "",
            "dna_alt_fraction": "",
            "decision": decision,
            "original_decision": "",
            "rescue_reason": "",
            "reference_id": "",
            "reference_gene": "",
            "reference_product": "",
            "reference_protein_id": "",
            "reference_inference": "",
            "protein_identity": "",
            "protein_coverage": "",
            "applied": False,
            "note": note,
        }
        evidence_rows.append(evidence_row)
        site["evidence_row"] = evidence_row
        site_records.append(site)
        if decision == "pending_dna":
            normal_dna_candidate_sites.append(site)

    ensure_dna_counts_for_sites(
        normal_dna_candidate_sites,
        dna_pileups=dna_pileups,
        record_name=record_name,
    )
    for site in normal_dna_candidate_sites:
        decision = classify_decision(
            rna_depth=site["rna_depth"],
            edited_reads=site["rna_edited_reads"],
            edit_fraction=site["rna_edit_fraction"],
            rna_mean_baseq=site["rna_mean_baseq"],
            rna_mean_mapq=site["rna_mean_mapq"],
            dna_depth=site["dna_depth"],
            dna_alt_fraction=site["dna_alt_fraction"],
            args=args,
        )
        site["evidence_row"]["decision"] = decision
        applied = decision == "accept"
        site["evidence_row"]["applied"] = applied
        if decision == "likely_genomic_variant":
            site["evidence_row"]["note"] = (
                "DNA/HiFi alternate fraction exceeds the RNA-editing "
                "threshold; review genomic sequence or variant."
            )
        elif decision == "insufficient_dna_coverage":
            site["evidence_row"]["note"] = (
                "RNA support passed normal candidate thresholds, but DNA/HiFi "
                "coverage was below the configured minimum."
            )
        else:
            site["evidence_row"]["note"] = ""

        summary["accepted_sites"] += int(decision == "accept")
        summary["applied_sites"] += int(applied)
        summary["likely_genomic_variant_sites"] += int(
            decision == "likely_genomic_variant"
        )
        if applied:
            edited_cds[site["cds_index"]] = "T"
            applied_indices.add(site["cds_index"])
            editing_effects.append(site["effect"])
            applied_decisions.append(decision)
            feature_decision["accepted_sites"].append(
                make_accepted_site(site, decision=decision)
            )

    if rna_depths:
        summary["mean_candidate_rna_depth"] = f"{statistics.fmean(rna_depths):.3f}"

    edited_coding_sequence = "".join(edited_cds)[codon_offset:]
    validation = validate_complete_cds_translation(edited_coding_sequence, table_id)
    if validation["valid"]:
        translation = translate_cds_for_genbank(edited_coding_sequence, table_id)
        add_exception = bool(applied_indices) and any(
            effect != "synonymous" for effect in editing_effects
        )
        updated_lines = update_cds_block_lines(
            block_lines,
            translation=translation,
            table_id=table_id,
            add_rna_editing_exception=add_exception,
            recalculate_terminal_notes=recalculate_terminal_notes,
            rna_editing_note=rna_editing_cds_note(editing_effects, applied_decisions),
            rna_editing_inferences=[RNA_EDITING_INFERENCE],
        )
        summary["translation_added"] = True
        summary["exception_added"] = add_exception and not any(
            value.lower() == "rna editing" for value in exceptions
        )
        summary["translation_status"] = "updated"
        feature_decision["translation_added"] = True
        feature_decision["exception_added"] = summary["exception_added"]
        feature_decision["translation_status"] = summary["translation_status"]
        return updated_lines, evidence_rows, summary, feature_decision

    ensure_dna_counts_for_sites(
        essential_rescue_dna_candidate_sites(
            validation["errors"],
            site_records,
            complete_coding_length // 3,
            args,
        ),
        dna_pileups=dna_pileups,
        record_name=record_name,
    )
    rescue_sites = essential_rescue_sites(
        validation["errors"],
        site_records,
        complete_coding_length // 3,
        args,
    )
    if rescue_sites:
        rescued_cds = list(edited_cds)
        for site, _reason in rescue_sites:
            rescued_cds[site["cds_index"]] = "T"
        rescued_coding_sequence = "".join(rescued_cds)[codon_offset:]
        rescued_validation = validate_complete_cds_translation(
            rescued_coding_sequence,
            table_id,
        )
        if rescued_validation["valid"]:
            for site, rescue_reason in rescue_sites:
                edited_cds[site["cds_index"]] = "T"
                applied_indices.add(site["cds_index"])
                editing_effects.append(site["effect"])
                applied_decisions.append("essential_rescue")
                evidence_row = site["evidence_row"]
                evidence_row["original_decision"] = evidence_row["decision"]
                evidence_row["decision"] = "essential_rescue"
                evidence_row["rescue_reason"] = rescue_reason
                evidence_row["applied"] = True
                evidence_row["note"] = (
                    "Applied only to rescue complete CDS validation."
                )
                feature_decision["accepted_sites"].append(
                    make_accepted_site(
                        site,
                        decision="essential_rescue",
                        rescue_reason=rescue_reason,
                    )
                )
                summary["applied_sites"] += 1
                summary["rescued_sites"] += 1

            translation = translate_cds_for_genbank(
                rescued_coding_sequence,
                table_id,
            )
            add_exception = bool(applied_indices) and any(
                effect != "synonymous" for effect in editing_effects
            )
            updated_lines = update_cds_block_lines(
                block_lines,
                translation=translation,
                table_id=table_id,
                add_rna_editing_exception=add_exception,
                recalculate_terminal_notes=recalculate_terminal_notes,
                rna_editing_note=rna_editing_cds_note(
                    editing_effects,
                    applied_decisions,
                ),
                rna_editing_inferences=[RNA_EDITING_INFERENCE],
            )
            summary["translation_added"] = True
            summary["exception_added"] = add_exception and not any(
                value.lower() == "rna editing" for value in exceptions
            )
            summary["translation_status"] = "updated_with_essential_rescue"
            feature_decision["translation_added"] = True
            feature_decision["exception_added"] = summary["exception_added"]
            feature_decision["translation_status"] = summary["translation_status"]
            return updated_lines, evidence_rows, summary, feature_decision

    reference_rescue_sites = []
    if references_by_gene is not None:
        summary["reference_inferred_rescue_attempted"] = True
        reference_rescue_sites = reference_inferred_candidate_sites(
            validation["errors"],
            site_records,
            complete_coding_length // 3,
        )
        summary["reference_inferred_candidate_sites"] = len(reference_rescue_sites)
    reference_rescue = None
    if references_by_gene is not None:
        reference_rescue = reference_inferred_rescue(
            validation_errors=validation["errors"],
            site_records=site_records,
            edited_cds=edited_cds,
            codon_offset=codon_offset,
            terminal_codon_index=complete_coding_length // 3,
            table_id=table_id,
            gene=gene,
            product=product,
            references_by_gene=references_by_gene,
            references_by_product=references_by_product,
            rescue_sites=reference_rescue_sites,
        )
    if reference_rescue is not None:
        reference = reference_rescue["reference"]
        metrics = reference_rescue["metrics"]
        reference_inference = reference_rescue["inference"]
        if not summary["product_filled"]:
            product_fill = reference_product_fill_from_reference(
                gene=gene,
                product_values=product_values,
                reference=reference,
                metrics=metrics,
            )
            if product_fill is not None:
                product = product_fill["product"]
                product_values = [product]
                block_lines = fill_cds_product_qualifier(block_lines, product)
                summary["product"] = product
                summary.update(product_fill_reference_summary(product_fill))
                feature_decision["product"] = product
                feature_decision["product_filled"] = True
                feature_decision["product_fill_reference"] = (
                    product_fill_reference_decision(product_fill)
                )
                for evidence_row in evidence_rows:
                    evidence_row["product"] = product
        for site, rescue_reason in reference_rescue["sites"]:
            edited_cds[site["cds_index"]] = "T"
            applied_indices.add(site["cds_index"])
            editing_effects.append(site["effect"])
            applied_decisions.append("reference_inferred")
            evidence_row = site["evidence_row"]
            was_likely_genomic_variant = (
                evidence_row["decision"] == "likely_genomic_variant"
            )
            evidence_row["original_decision"] = evidence_row["decision"]
            evidence_row["decision"] = "reference_inferred"
            evidence_row["rescue_reason"] = rescue_reason
            evidence_row["reference_id"] = reference.get("record_id", "")
            evidence_row["reference_gene"] = reference.get("gene", "")
            evidence_row["reference_product"] = reference.get("product", "")
            evidence_row["reference_protein_id"] = reference.get("protein_id", "")
            evidence_row["reference_inference"] = reference_inference
            evidence_row["protein_identity"] = f"{metrics['protein_identity']:.6f}"
            evidence_row["protein_coverage"] = f"{metrics['protein_coverage']:.6f}"
            evidence_row["applied"] = True
            evidence_row["note"] = (
                "Applied only to rescue complete CDS validation using "
                "reference protein evidence."
            )
            feature_decision["accepted_sites"].append(
                make_accepted_site(
                    site,
                    decision="reference_inferred",
                    rescue_reason=rescue_reason,
                    reference=reference,
                    metrics=metrics,
                )
            )
            summary["applied_sites"] += 1
            summary["reference_inferred_sites"] += 1
            summary["likely_genomic_variant_sites"] -= int(was_likely_genomic_variant)

        add_exception = bool(applied_indices) and any(
            effect != "synonymous" for effect in editing_effects
        )
        inferences = []
        if any(
            decision in {"accept", "essential_rescue"}
            for decision in applied_decisions
        ):
            inferences.append(RNA_EDITING_INFERENCE)
        inferences.append(reference_inference)
        updated_lines = update_cds_block_lines(
            block_lines,
            translation=reference_rescue["translation"],
            table_id=table_id,
            add_rna_editing_exception=add_exception,
            recalculate_terminal_notes=recalculate_terminal_notes,
            rna_editing_note=rna_editing_cds_note(
                editing_effects,
                applied_decisions,
            ),
            rna_editing_inferences=inferences,
        )
        summary["translation_added"] = True
        summary["exception_added"] = add_exception and not any(
            value.lower() == "rna editing" for value in exceptions
        )
        summary["translation_status"] = "updated_with_reference_inference"
        feature_decision["translation_added"] = True
        feature_decision["exception_added"] = summary["exception_added"]
        feature_decision["translation_status"] = summary["translation_status"]
        return updated_lines, evidence_rows, summary, feature_decision

    undetermined_terminal_validation = (
        validate_cds_translation_allowing_undetermined_terminal(
            edited_coding_sequence,
            table_id,
        )
    )
    if undetermined_terminal_validation["valid"]:
        start_undetermined = undetermined_terminal_validation["start_undetermined"]
        stop_undetermined = undetermined_terminal_validation["stop_undetermined"]
        translation = translate_cds_for_genbank(edited_coding_sequence, table_id)
        updated_location = add_partial_markers_to_location(
            block.location,
            five_prime=start_undetermined,
            three_prime=stop_undetermined,
        )
        extra_notes = []
        if start_undetermined:
            extra_notes.append(START_CODON_UNDETERMINED_NOTE)
        if stop_undetermined:
            extra_notes.append(STOP_CODON_UNDETERMINED_NOTE)
        add_exception = bool(applied_indices) and any(
            effect != "synonymous" for effect in editing_effects
        )
        updated_lines = update_cds_block_lines(
            block_lines,
            translation=translation,
            table_id=table_id,
            add_rna_editing_exception=add_exception,
            location=updated_location,
            recalculate_terminal_notes=recalculate_terminal_notes,
            rna_editing_note=rna_editing_cds_note(
                editing_effects,
                applied_decisions,
            ),
            rna_editing_inferences=[RNA_EDITING_INFERENCE],
            extra_notes=extra_notes,
        )
        summary["location"] = updated_location
        summary["translation_added"] = True
        summary["exception_added"] = add_exception and not any(
            value.lower() == "rna editing" for value in exceptions
        )
        summary["translation_status"] = undetermined_terminal_translation_status(
            start_undetermined,
            stop_undetermined,
        )
        feature_decision["location"] = updated_location
        feature_decision["translation_added"] = True
        feature_decision["exception_added"] = summary["exception_added"]
        feature_decision["translation_status"] = summary["translation_status"]
        return updated_lines, evidence_rows, summary, feature_decision

    summary["translation_status"] = "validation_failed"
    summary["translation_errors"] = ";".join(validation["errors"])
    feature_decision["translation_status"] = summary["translation_status"]
    feature_decision["translation_errors"] = validation["errors"]
    return block_lines, evidence_rows, summary, feature_decision


def curate_records(args, rna_bams, dna_bams):
    lines = args.input_gbk.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    curated_records = []
    evidence_rows = []
    summary_rows = []
    feature_decisions = []
    rna_pileups = PileupBatchCounter(
        rna_bams,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
    )
    dna_pileups = PileupBatchCounter(
        dna_bams,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
    )
    references_by_gene = None
    references_by_product = None
    if args.reference_dir is not None:
        if not args.reference_dir.is_dir():
            raise RuntimeError(
                f"Reference directory does not exist: {args.reference_dir}"
            )
        references_by_gene = reference_cds_by_gene(
            args.reference_dir,
            args.transl_table,
        )
        references_by_product = reference_cds_by_product(references_by_gene)

    for record_index, record in enumerate(records, start=1):
        metadata = genbank_record_locus_metadata(record)
        record_name = (
            metadata["name"] if metadata is not None else f"record_{record_index}"
        )
        record_sequence = genbank_record_origin_sequence(record)
        features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
        if features_index is None or origin_index is None:
            curated_records.append(record)
            continue

        rna_pileups.preload(
            record_name,
            candidate_positions_for_record(blocks, record_sequence),
            full_span=True,
        )
        curated_blocks = []
        (
            existing_site_features,
            existing_site_feature_lines,
        ) = existing_rna_editing_site_feature_registry(blocks)
        for block in blocks:
            if block.key != "CDS":
                curated_blocks.append(block.lines)
                continue
            (
                updated_lines,
                block_evidence_rows,
                summary,
                feature_decision,
            ) = evaluate_cds_feature(
                record_name=record_name,
                record_sequence=record_sequence,
                block=block,
                args=args,
                rna_pileups=rna_pileups,
                dna_pileups=dna_pileups,
                references_by_gene=references_by_gene,
                references_by_product=references_by_product,
            )
            curated_blocks.append(updated_lines)
            curated_blocks.extend(
                new_rna_editing_site_feature_blocks(
                    feature_decision,
                    existing_site_features,
                    existing_site_feature_lines,
                )
            )
            evidence_rows.extend(block_evidence_rows)
            summary_rows.append(summary)
            feature_decisions.append(feature_decision)

        curated_feature_lines = [
            line for block_lines in curated_blocks for line in block_lines
        ]
        curated_records.append(
            [
                *record[: features_index + 1],
                *curated_feature_lines,
                *record[origin_index:],
            ]
        )

    curated_lines = [line for record in curated_records for line in record]
    return curated_lines, evidence_rows, summary_rows, feature_decisions


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def count_summary(summary_rows, key):
    return sum(int(row.get(key, 0) or 0) for row in summary_rows)


def count_true(summary_rows, key):
    return sum(1 for row in summary_rows if bool(row.get(key)))


def format_path_for_markdown(path: Path):
    return f"`{path.as_posix()}`"


def format_applied_site_count_breakdown(
    accepted_site_count,
    applied_site_count,
    rescued_site_count,
    reference_inferred_site_count,
):
    parts = [f"{accepted_site_count} normal-threshold accepted site(s)"]
    if rescued_site_count:
        parts.append(f"{rescued_site_count} CDS-essential rescue-applied site(s)")
    if reference_inferred_site_count:
        parts.append(
            f"{reference_inferred_site_count} reference-inferred rescue-applied "
            "site(s)"
        )
    parts.append(f"{applied_site_count} total applied site(s)")
    return ", ".join(parts)


def format_rna_editing_post_curation_section(
    args,
    evidence_rows,
    summary_rows,
    cds_qualifier_order=None,
):
    cds_count = len(summary_rows)
    candidate_site_count = count_summary(summary_rows, "candidate_c_sites")
    accepted_site_count = count_summary(summary_rows, "accepted_sites")
    applied_site_count = count_summary(summary_rows, "applied_sites")
    rescued_site_count = count_summary(summary_rows, "rescued_sites")
    reference_inferred_site_count = count_summary(
        summary_rows,
        "reference_inferred_sites",
    )
    reference_inferred_attempted_cds_count = count_true(
        summary_rows,
        "reference_inferred_rescue_attempted",
    )
    reference_inferred_candidate_site_count = count_summary(
        summary_rows,
        "reference_inferred_candidate_sites",
    )
    reference_inferred_fixed_cds_count = sum(
        1
        for row in summary_rows
        if row.get("translation_status") == "updated_with_reference_inference"
    )
    likely_genomic_variant_site_count = count_summary(
        summary_rows,
        "likely_genomic_variant_sites",
    )
    translated_cds_count = count_true(summary_rows, "translation_added")
    exception_cds_count = count_true(summary_rows, "exception_added")
    product_filled_cds_count = count_true(summary_rows, "product_filled")
    _ = evidence_rows
    validation_failed_count = sum(
        1
        for row in summary_rows
        if row.get("translation_status") == "validation_failed"
    )
    sidecars = ", ".join(
        format_path_for_markdown(path)
        for path in (args.evidence_tsv, args.decisions_json, args.summary_tsv)
    )
    rna_bams = ", ".join(format_path_for_markdown(path) for path in args.rna_bam)
    dna_bam = format_path_for_markdown(args.dna_bam) if args.dna_bam else "none"
    applied_site_breakdown = format_applied_site_count_breakdown(
        accepted_site_count,
        applied_site_count,
        rescued_site_count,
        reference_inferred_site_count,
    )
    translation_rescue_context = (
        "normal RNA-editing calls, CDS-essential rescue, and "
        "reference-inferred rescue"
        if args.reference_dir is not None
        else "normal RNA-editing calls and CDS-essential rescue"
    )

    lines = [
        "## RNA editing evidence curation",
        "",
        (
            "- RNA editing curation: scanned all annotated CDS for plant "
            "organelle C-to-U editing candidates and regenerated CDS "
            "`/translation` qualifiers from edited CDS sequences where evidence "
            "passed the accept thresholds."
        ),
        (
            f"- RNA editing evidence: {candidate_site_count} candidate C-to-U "
            f"site(s) across {cds_count} CDS feature(s); {applied_site_breakdown}."
        ),
        (
            f"- RNA editing GenBank updates: added or refreshed `/translation` "
            f"for {translated_cds_count} CDS feature(s); added "
            f"`/exception=\"RNA editing\"` for {exception_cds_count} CDS "
            "feature(s), with DDBJ-oriented CDS `/note` and `/inference` "
            "qualifiers where RNA editing changes the translated product. "
            "Applied edit sites are also represented as `misc_feature` "
            f"annotations with `/note=\"{RNA_EDITING_SITE_NOTE}\"`; "
            "reference-inferred sites include an additional reference CDS "
            "evidence note."
        ),
        (
            "- RNA editing thresholds: accept requires "
            f"RNA depth >= {args.min_rna_depth}, edited reads >= "
            f"{args.min_edited_reads}, edit fraction >= "
            f"{args.min_edit_fraction:g}, base quality >= "
            f"{args.min_base_quality}, mapping quality >= "
            f"{args.min_mapping_quality}, DNA/HiFi depth >= "
            f"{args.min_dna_depth}, and DNA alternate fraction <= "
            f"{args.max_dna_alt_fraction:g}."
        ),
        (
            "- RNA editing non-accepted calls: retained in the evidence "
            "tables but not applied to the GenBank output, except for "
            "CDS-essential rescue sites that restore a valid start codon, "
            "terminal stop codon, or internal premature-stop rescue, and "
            "reference-inferred sites that rescue otherwise invalid CDS "
            "translation. DNA/HiFi pileup is evaluated only for sites with "
            "normal RNA-editing candidate support or CDS-essential rescue "
            "relevance; reference-inferred rescue does not require HiFi/DNA "
            "support."
        ),
        (
            "- RNA editing CDS-essential rescue thresholds: RNA depth >= "
            f"{args.essential_rescue_min_rna_depth}, edited reads >= "
            f"{args.essential_rescue_min_edited_reads}, edit fraction >= "
            f"{args.essential_rescue_min_edit_fraction:g}, DNA/HiFi depth >= "
            f"{args.min_dna_depth}, and DNA alternate fraction <= "
            f"{args.essential_rescue_max_dna_alt_fraction:g}."
        ),
        f"- RNA editing RNA-seq BAM input(s): {rna_bams}.",
        f"- RNA editing DNA/HiFi BAM input: {dna_bam}.",
        f"- RNA editing sidecars: {sidecars}.",
    ]
    if likely_genomic_variant_site_count:
        lines.append(
            "- RNA editing genome-sequence review note: "
            f"{likely_genomic_variant_site_count} candidate site(s) had "
            "DNA/HiFi alternate support above the configured threshold and "
            "were classified as `likely_genomic_variant`; these were retained "
            "in the evidence table and not applied as RNA editing."
        )
    if cds_qualifier_order:
        target_feature_count = cds_qualifier_order[
            "cds_qualifier_order_target_feature_count"
        ]
        changed_feature_count = cds_qualifier_order[
            "cds_qualifier_order_changed_feature_count"
        ]
        changed_record_count = cds_qualifier_order[
            "cds_qualifier_order_changed_record_count"
        ]
        record_count = cds_qualifier_order["cds_qualifier_order_record_count"]
        if cds_qualifier_order["cds_qualifier_order_changed"]:
            lines.append(
                "- RNA editing GenBank CDS qualifier order: normalized "
                f"{changed_feature_count}/{target_feature_count} CDS feature(s) "
                f"across {changed_record_count}/{record_count} record(s) using "
                "a stable qualifier-block sort; repeated and unranked "
                "qualifiers kept their original relative order."
            )
        else:
            lines.append(
                "- RNA editing GenBank CDS qualifier order: already normalized "
                f"across {target_feature_count} CDS feature(s) in "
                f"{record_count} record(s)."
            )
    if args.reference_dir is not None:
        reference_dir = format_path_for_markdown(args.reference_dir)
        if product_filled_cds_count:
            lines.append(
                "- RNA editing reference product fill: filled "
                f"{product_filled_cds_count} empty CDS `/product` qualifier(s) "
                "from same-gene reference CDS/protein matches with "
                "identity/coverage >= "
                f"{REFERENCE_PRODUCT_FILL_MIN_PROTEIN_IDENTITY:g}/"
                f"{REFERENCE_PRODUCT_FILL_MIN_PROTEIN_COVERAGE:g}; existing "
                "non-empty `/product` qualifiers were left unchanged."
            )
        lines.append(
            "- RNA editing reference-inferred rescue: "
            f"enabled with reference directory {reference_dir}; evaluated "
            f"{reference_inferred_attempted_cds_count} CDS feature(s) that "
            "still failed complete-CDS validation after normal RNA-editing "
            f"calls and CDS-essential rescue; tested "
            f"{reference_inferred_candidate_site_count} validation-targeted "
            f"candidate site(s); fixed {reference_inferred_fixed_cds_count} "
            f"CDS feature(s) by applying {reference_inferred_site_count} "
            "reference-inferred site(s). Candidate sites are selected only "
            "when exactly one not-yet-applied C-to-U site directly addresses "
            "a remaining validation error (start gain at the first codon, "
            "terminal stop gain, or internal premature-stop rescue). A fix is "
            "applied only if the edited CDS passes complete-CDS translation validation "
            "and the translated protein matches the best same-gene/product "
            "reference with identity/coverage >= "
            f"{REFERENCE_INFERRED_MIN_PROTEIN_IDENTITY:g}/"
            f"{REFERENCE_INFERRED_MIN_PROTEIN_COVERAGE:g}."
        )
    if validation_failed_count:
        lines.append(
            "- RNA editing translation note: after "
            f"{translation_rescue_context}, {validation_failed_count} CDS "
            "feature(s) were still not updated because the edited CDS did "
            "not pass complete-CDS translation validation."
        )
    return "\n".join(lines)


def insert_post_curation_section(text, section):
    marker = "\n# Manual post-curation"
    section = section.strip()
    if marker in text:
        return text.replace(marker, f"\n\n{section}\n{marker}", 1)
    if text.strip():
        return f"{text.rstrip()}\n\n{section}\n"
    return (
        "# Automatic post-curation\n\n"
        f"{section}\n\n"
        "# Manual post-curation\n"
    )


def write_post_curation(args, evidence_rows, summary_rows, cds_qualifier_order=None):
    if args.output_post_curation is None:
        return
    if args.input_post_curation is not None and args.input_post_curation.exists():
        text = args.input_post_curation.read_text()
    else:
        text = ""
    section = format_rna_editing_post_curation_section(
        args,
        evidence_rows,
        summary_rows,
        cds_qualifier_order=cds_qualifier_order,
    )
    args.output_post_curation.parent.mkdir(parents=True, exist_ok=True)
    args.output_post_curation.write_text(insert_post_curation_section(text, section))


def main():
    args = parse_args()
    rna_bams = open_bams(args.rna_bam)
    dna_bams = open_bams([args.dna_bam]) if args.dna_bam else []
    try:
        curated_lines, evidence_rows, summary_rows, feature_decisions = curate_records(
            args,
            rna_bams,
            dna_bams,
        )
    finally:
        close_bams(rna_bams)
        close_bams(dna_bams)

    args.output_gbk.parent.mkdir(parents=True, exist_ok=True)
    args.output_gbk.write_text("".join(curated_lines))
    cds_qualifier_order = normalize_genbank_cds_qualifier_order(args.output_gbk)
    write_tsv(args.evidence_tsv, evidence_rows, EVIDENCE_FIELDS)
    write_tsv(args.summary_tsv, summary_rows, SUMMARY_FIELDS)
    write_post_curation(
        args,
        evidence_rows,
        summary_rows,
        cds_qualifier_order=cds_qualifier_order,
    )

    decisions = {
        "input_gbk": str(args.input_gbk),
        "output_gbk": str(args.output_gbk),
        "organelle": args.organelle,
        "tool": args.tool,
        "transl_table_default": str(args.transl_table),
        "reference_dir": str(args.reference_dir) if args.reference_dir else None,
        "rna_bams": [str(path) for path in args.rna_bam],
        "dna_bam": str(args.dna_bam) if args.dna_bam else None,
        "thresholds": {
            "min_rna_depth": args.min_rna_depth,
            "min_edited_reads": args.min_edited_reads,
            "min_edit_fraction": args.min_edit_fraction,
            "essential_rescue_min_rna_depth": args.essential_rescue_min_rna_depth,
            "essential_rescue_min_edited_reads": args.essential_rescue_min_edited_reads,
            "essential_rescue_min_edit_fraction": args.essential_rescue_min_edit_fraction,
            "essential_rescue_max_dna_alt_fraction": args.essential_rescue_max_dna_alt_fraction,
            "min_base_quality": args.min_base_quality,
            "min_mapping_quality": args.min_mapping_quality,
            "min_dna_depth": args.min_dna_depth,
            "max_dna_alt_fraction": args.max_dna_alt_fraction,
            "reference_inferred_min_protein_identity": (
                REFERENCE_INFERRED_MIN_PROTEIN_IDENTITY
            ),
            "reference_inferred_min_protein_coverage": (
                REFERENCE_INFERRED_MIN_PROTEIN_COVERAGE
            ),
            "reference_product_fill_min_protein_identity": (
                REFERENCE_PRODUCT_FILL_MIN_PROTEIN_IDENTITY
            ),
            "reference_product_fill_min_protein_coverage": (
                REFERENCE_PRODUCT_FILL_MIN_PROTEIN_COVERAGE
            ),
        },
        "cds_qualifier_order": cds_qualifier_order,
        "feature_decisions": feature_decisions,
    }
    args.decisions_json.parent.mkdir(parents=True, exist_ok=True)
    args.decisions_json.write_text(
        json.dumps(decisions, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
