#!/usr/bin/env python3
import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

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
    genbank_feature_has_qualifier,
    genbank_qualifier_name,
    genbank_qualifier_values,
    genbank_record_locus_metadata,
    genbank_record_origin_sequence,
    parse_genbank_feature_blocks,
    split_genbank_location_arguments,
    split_genbank_records,
    translate_complete_codons,
    unwrap_genbank_location_function,
    validate_complete_cds_translation,
)


BASES = ("A", "C", "G", "T")
RNA_EDITING_EXCEPTION = "RNA editing"
RNA_EDITING_INFERENCE = "EXISTENCE:similar to RNA sequence, mRNA (same species)"
RNA_EDITING_SITE_NOTE = "C to U RNA editing"
RNA_EDITING_TRANSLATION_NOTE = "RNA editing is required for translation"
RNA_EDITING_START_GAIN_NOTE = "initiation codon is created by RNA editing"
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
    "applied",
    "note",
]
SUMMARY_FIELDS = [
    "record",
    "organelle",
    "tool",
    "gene",
    "product",
    "feature_index",
    "location",
    "candidate_c_sites",
    "high_sites",
    "moderate_sites",
    "applied_sites",
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
    parser.add_argument("--min-rna-depth", type=int, default=20)
    parser.add_argument("--min-edited-reads", type=int, default=5)
    parser.add_argument("--min-edit-fraction", type=float, default=0.20)
    parser.add_argument("--moderate-min-rna-depth", type=int, default=10)
    parser.add_argument("--moderate-min-edited-reads", type=int, default=3)
    parser.add_argument("--moderate-min-edit-fraction", type=float, default=0.10)
    parser.add_argument("--min-base-quality", type=int, default=30)
    parser.add_argument("--min-mapping-quality", type=int, default=30)
    parser.add_argument("--min-dna-depth", type=int, default=10)
    parser.add_argument("--max-dna-alt-fraction", type=float, default=0.05)
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
    high = (
        rna_depth >= args.min_rna_depth
        and edited_reads >= args.min_edited_reads
        and edit_fraction >= args.min_edit_fraction
        and rna_mean_baseq >= args.min_base_quality
        and rna_mean_mapq >= args.min_mapping_quality
        and dna_supported
    )
    if high:
        return "high"

    moderate = (
        rna_depth >= args.moderate_min_rna_depth
        and edited_reads >= args.moderate_min_edited_reads
        and edit_fraction >= args.moderate_min_edit_fraction
        and dna_supported
    )
    if moderate:
        return "moderate"
    if rna_depth == 0:
        return "insufficient_coverage"
    if dna_depth < args.min_dna_depth:
        return "insufficient_dna_coverage"
    if dna_alt_fraction > args.max_dna_alt_fraction:
        return "likely_genomic_variant"
    return "reject"


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


def update_cds_block_lines(
    block_lines,
    *,
    translation,
    table_id,
    add_rna_editing_exception,
    rna_editing_note=None,
):
    cleaned = remove_qualifiers(block_lines, {"translation"})
    existing_exceptions = genbank_qualifier_values(cleaned, "exception")
    has_rna_editing_exception = any(
        value.lower() == RNA_EDITING_EXCEPTION.lower() for value in existing_exceptions
    )
    if add_rna_editing_exception and not has_rna_editing_exception:
        cleaned.extend(format_quoted_qualifier("exception", RNA_EDITING_EXCEPTION))
    if add_rna_editing_exception:
        existing_notes = {
            value.lower() for value in genbank_qualifier_values(cleaned, "note")
        }
        if rna_editing_note and rna_editing_note.lower() not in existing_notes:
            cleaned.extend(format_quoted_qualifier("note", rna_editing_note))
        existing_inferences = {
            value.lower()
            for value in genbank_qualifier_values(cleaned, "inference")
        }
        if RNA_EDITING_INFERENCE.lower() not in existing_inferences:
            cleaned.extend(format_quoted_qualifier("inference", RNA_EDITING_INFERENCE))
    if table_id != "1" and not any(genbank_qualifier_name(line) == "transl_table" for line in cleaned):
        cleaned.extend(format_simple_qualifier("transl_table", f"\"{table_id}\""))
    cleaned.extend(format_quoted_qualifier("translation", translation))
    return cleaned


def rna_editing_cds_note(editing_effects):
    if "start_gain" in editing_effects:
        return RNA_EDITING_START_GAIN_NOTE
    return RNA_EDITING_TRANSLATION_NOTE


def rna_editing_site_location(site):
    position = site["genomic_pos"]
    if site.get("genomic_base") == "G":
        return f"complement({position})"
    return str(position)


def rna_editing_site_feature_lines(site):
    location = rna_editing_site_location(site)
    return [
        f"{FEATURE_INDENT}{'misc_feature':<16}{location}\n",
        *format_quoted_qualifier("note", RNA_EDITING_SITE_NOTE),
    ]


def existing_rna_editing_site_features(blocks):
    existing = set()
    for block in blocks:
        if block.key != "misc_feature":
            continue
        notes = {
            value.lower() for value in genbank_qualifier_values(block.lines, "note")
        }
        if RNA_EDITING_SITE_NOTE.lower() in notes:
            existing.add((block.location, RNA_EDITING_SITE_NOTE.lower()))
    return existing


def new_rna_editing_site_feature_blocks(feature_decision, existing_site_features):
    blocks = []
    for site in feature_decision["accepted_sites"]:
        location = rna_editing_site_location(site)
        key = (location, RNA_EDITING_SITE_NOTE.lower())
        if key in existing_site_features:
            continue
        existing_site_features.add(key)
        blocks.append(rna_editing_site_feature_lines(site))
    return blocks


def feature_candidate_sources(tool, exceptions):
    sources = ["all_CDS_scan"]
    if any(value.lower() == "rna editing" for value in exceptions):
        sources.append(f"{tool}_existing_RNA_editing_exception")
    return ";".join(sources)


def empty_counts():
    return {base: 0 for base in BASES}


def evaluate_cds_feature(
    *,
    record_name,
    record_sequence,
    block,
    args,
    rna_bams,
    dna_bams,
):
    table_id = str(
        first_genbank_qualifier_value(block.lines, "transl_table", args.transl_table)
    )
    if table_id not in GENETIC_CODES:
        table_id = str(args.transl_table)
    gene = first_genbank_qualifier_value(block.lines, "gene", "?")
    product = first_genbank_qualifier_value(block.lines, "product", "")
    exceptions = genbank_qualifier_values(block.lines, "exception")
    candidate_source = feature_candidate_sources(args.tool, exceptions)
    is_pseudo = genbank_feature_has_qualifier(
        block.lines,
        "pseudo",
    ) or genbank_feature_has_qualifier(block.lines, "pseudogene")
    is_partial = "<" in block.location or ">" in block.location
    strand = infer_location_strand(block.location)

    summary = {
        "record": record_name,
        "organelle": args.organelle,
        "tool": args.tool,
        "gene": gene,
        "product": product,
        "feature_index": block.feature_index,
        "location": block.location,
        "candidate_c_sites": 0,
        "high_sites": 0,
        "moderate_sites": 0,
        "applied_sites": 0,
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
        "location": block.location,
        "candidate_source": candidate_source,
        "accepted_sites": [],
        "translation_added": False,
        "exception_added": False,
        "translation_status": "not_evaluated",
        "translation_errors": [],
    }

    evidence_rows = []
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

    edited_cds = list(cds_sequence)
    rna_depths = []
    applied_indices = set()
    editing_effects = []
    complete_coding_length = len(cds_sequence[codon_offset:])
    complete_coding_length -= complete_coding_length % 3
    coding_end = codon_offset + complete_coding_length

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
        rna_counts = pileup_base_counts(
            rna_bams,
            record_name,
            genomic_pos - 1,
            min_base_quality=args.min_base_quality,
            min_mapping_quality=args.min_mapping_quality,
        )
        if dna_bams:
            dna_counts = pileup_base_counts(
                dna_bams,
                record_name,
                genomic_pos - 1,
                min_base_quality=args.min_base_quality,
                min_mapping_quality=args.min_mapping_quality,
            )
        else:
            dna_counts = {
                "counts": empty_counts(),
                "depth": 0,
                "mean_baseq": 0.0,
                "mean_mapq": 0.0,
                "sample_depths": [],
            }
        rna_depth = rna_counts["depth"]
        rna_edited_reads = rna_counts["counts"][edited_aligned_base]
        rna_edit_fraction = (
            rna_edited_reads / rna_depth if rna_depth > 0 else 0.0
        )
        dna_depth = dna_counts["depth"]
        dna_alt_reads = dna_counts["counts"][edited_aligned_base]
        dna_alt_fraction = dna_alt_reads / dna_depth if dna_depth > 0 else 1.0
        decision = classify_decision(
            rna_depth=rna_depth,
            edited_reads=rna_edited_reads,
            edit_fraction=rna_edit_fraction,
            rna_mean_baseq=rna_counts["mean_baseq"],
            rna_mean_mapq=rna_counts["mean_mapq"],
            dna_depth=dna_depth,
            dna_alt_fraction=dna_alt_fraction,
            args=args,
        )
        applied = decision == "high"
        if applied:
            edited_cds[cds_index] = "T"
            applied_indices.add(cds_index)
            editing_effects.append(effect)
            feature_decision["accepted_sites"].append(
                {
                    "genomic_pos": genomic_pos,
                    "cds_pos": cds_index + 1,
                    "codon_index": codon_index,
                    "codon_pos": codon_pos,
                    "genomic_base": genomic_base,
                    "edited_base": edited_aligned_base,
                    "effect": effect,
                    "decision": decision,
                    "rna_depth": rna_depth,
                    "rna_edited_reads": rna_edited_reads,
                    "rna_edit_fraction": rna_edit_fraction,
                    "dna_depth": dna_depth,
                    "dna_alt_fraction": dna_alt_fraction,
                }
            )

        rna_depths.append(rna_depth)
        summary["candidate_c_sites"] += 1
        summary["high_sites"] += int(decision == "high")
        summary["moderate_sites"] += int(decision == "moderate")
        summary["applied_sites"] += int(applied)
        summary["candidate_sites_with_min_depth"] += int(
            rna_depth >= args.min_rna_depth
        )
        evidence_rows.append(
            {
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
                "dna_depth": dna_depth,
                "dna_alt_reads": dna_alt_reads,
                "dna_alt_fraction": f"{dna_alt_fraction:.6f}",
                "decision": decision,
                "applied": applied,
                "note": "",
            }
        )

    if rna_depths:
        summary["mean_candidate_rna_depth"] = f"{statistics.fmean(rna_depths):.3f}"

    edited_coding_sequence = "".join(edited_cds)[codon_offset:]
    if is_partial:
        summary["translation_status"] = "skipped_partial"
        feature_decision["translation_status"] = summary["translation_status"]
        return block.lines, evidence_rows, summary, feature_decision

    validation = validate_complete_cds_translation(edited_coding_sequence, table_id)
    if validation["valid"]:
        translation = translate_cds_for_genbank(edited_coding_sequence, table_id)
        add_exception = bool(applied_indices) and any(
            effect != "synonymous" for effect in editing_effects
        )
        updated_lines = update_cds_block_lines(
            block.lines,
            translation=translation,
            table_id=table_id,
            add_rna_editing_exception=add_exception,
            rna_editing_note=rna_editing_cds_note(editing_effects),
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

    summary["translation_status"] = "validation_failed"
    summary["translation_errors"] = ";".join(validation["errors"])
    feature_decision["translation_status"] = summary["translation_status"]
    feature_decision["translation_errors"] = validation["errors"]
    return block.lines, evidence_rows, summary, feature_decision


def curate_records(args, rna_bams, dna_bams):
    lines = args.input_gbk.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    curated_records = []
    evidence_rows = []
    summary_rows = []
    feature_decisions = []

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

        curated_blocks = []
        existing_site_features = existing_rna_editing_site_features(blocks)
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
                rna_bams=rna_bams,
                dna_bams=dna_bams,
            )
            curated_blocks.append(updated_lines)
            curated_blocks.extend(
                new_rna_editing_site_feature_blocks(
                    feature_decision,
                    existing_site_features,
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


def format_rna_editing_post_curation_section(args, evidence_rows, summary_rows):
    cds_count = len(summary_rows)
    candidate_site_count = count_summary(summary_rows, "candidate_c_sites")
    high_site_count = count_summary(summary_rows, "high_sites")
    moderate_site_count = count_summary(summary_rows, "moderate_sites")
    applied_site_count = count_summary(summary_rows, "applied_sites")
    translated_cds_count = count_true(summary_rows, "translation_added")
    exception_cds_count = count_true(summary_rows, "exception_added")
    insufficient_coverage_count = sum(
        1
        for row in evidence_rows
        if row.get("decision") == "insufficient_coverage"
    )
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

    lines = [
        "## RNA editing evidence curation",
        "",
        (
            "- RNA editing curation: scanned all annotated CDS for plant "
            "organelle C-to-U editing candidates and regenerated CDS "
            "`/translation` qualifiers from edited CDS sequences where evidence "
            "passed the high-confidence thresholds."
        ),
        (
            f"- RNA editing evidence: {candidate_site_count} candidate C-to-U "
            f"site(s) across {cds_count} CDS feature(s); {high_site_count} "
            f"high-confidence site(s), {moderate_site_count} moderate site(s), "
            f"{applied_site_count} site(s) applied."
        ),
        (
            f"- RNA editing GenBank updates: added or refreshed `/translation` "
            f"for {translated_cds_count} CDS feature(s); added "
            f"`/exception=\"RNA editing\"` for {exception_cds_count} CDS "
            "feature(s), with DDBJ-oriented CDS `/note` and `/inference` "
            "qualifiers where RNA editing changes the translated product. "
            "Applied edit sites are also represented as `misc_feature` "
            f"annotations with `/note=\"{RNA_EDITING_SITE_NOTE}\"`."
        ),
        (
            "- RNA editing thresholds: high confidence requires "
            f"RNA depth >= {args.min_rna_depth}, edited reads >= "
            f"{args.min_edited_reads}, edit fraction >= "
            f"{args.min_edit_fraction:g}, base quality >= "
            f"{args.min_base_quality}, mapping quality >= "
            f"{args.min_mapping_quality}, DNA/HiFi depth >= "
            f"{args.min_dna_depth}, and DNA alternate fraction <= "
            f"{args.max_dna_alt_fraction:g}."
        ),
        (
            "- RNA editing moderate calls: retained in the evidence tables but "
            "not applied to the GenBank output."
        ),
        f"- RNA editing RNA-seq BAM input(s): {rna_bams}.",
        f"- RNA editing DNA/HiFi BAM input: {dna_bam}.",
        f"- RNA editing sidecars: {sidecars}.",
    ]
    if insufficient_coverage_count:
        lines.append(
            f"- RNA editing coverage note: {insufficient_coverage_count} "
            "candidate site(s) had no passing RNA-seq pileup coverage."
        )
    if validation_failed_count:
        lines.append(
            f"- RNA editing translation note: {validation_failed_count} CDS "
            "feature(s) were not updated because the edited CDS did not pass "
            "complete-CDS translation validation."
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


def write_post_curation(args, evidence_rows, summary_rows):
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
    write_tsv(args.evidence_tsv, evidence_rows, EVIDENCE_FIELDS)
    write_tsv(args.summary_tsv, summary_rows, SUMMARY_FIELDS)
    write_post_curation(args, evidence_rows, summary_rows)

    decisions = {
        "input_gbk": str(args.input_gbk),
        "output_gbk": str(args.output_gbk),
        "organelle": args.organelle,
        "tool": args.tool,
        "transl_table_default": str(args.transl_table),
        "rna_bams": [str(path) for path in args.rna_bam],
        "dna_bam": str(args.dna_bam) if args.dna_bam else None,
        "thresholds": {
            "min_rna_depth": args.min_rna_depth,
            "min_edited_reads": args.min_edited_reads,
            "min_edit_fraction": args.min_edit_fraction,
            "moderate_min_rna_depth": args.moderate_min_rna_depth,
            "moderate_min_edited_reads": args.moderate_min_edited_reads,
            "moderate_min_edit_fraction": args.moderate_min_edit_fraction,
            "min_base_quality": args.min_base_quality,
            "min_mapping_quality": args.min_mapping_quality,
            "min_dna_depth": args.min_dna_depth,
            "max_dna_alt_fraction": args.max_dna_alt_fraction,
        },
        "feature_decisions": feature_decisions,
    }
    args.decisions_json.parent.mkdir(parents=True, exist_ok=True)
    args.decisions_json.write_text(
        json.dumps(decisions, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
