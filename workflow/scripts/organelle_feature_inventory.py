#!/usr/bin/env python3
"""Read-only gene presence, copy-number, and duplicate-annotation review."""

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from reference_cds_qc import (
    REFERENCE_EXTENSIONS,
    SeqIO,
    feature_location_string,
    normalized_genbank_text_for_biopython,
    record_locus_id,
    require_biopython,
)


PRODUCT_TYPES = {"CDS", "tRNA", "rRNA"}
DUPLICATE_TYPES = PRODUCT_TYPES | {"gene"}
INVENTORY_FIELDS = [
    "gene", "target_feature_types", "target_copy_count", "target_feature_count",
    "target_pseudogene_count", "target_partial_count", "reference_count",
    "reference_present_count", "reference_copy_min", "reference_copy_max",
    "reference_copy_counts", "target_locations", "reference_comparison",
    "review_reason", "manual_review_candidate",
]
DUPLICATE_FIELDS = [
    "record_id", "gene", "feature_type", "location", "feature_indices",
    "locus_tags", "feature_count",
]


@dataclass(frozen=True)
class Feature:
    record_id: str
    record_index: int
    feature_index: int
    gene: str
    feature_type: str
    location: str
    location_key: tuple
    intervals: tuple
    locus_tag: str
    pseudogene: bool
    partial: bool

    @property
    def gene_key(self):
        return self.gene.casefold()

    @property
    def copy_key(self):
        # Coordinates on separate contigs are separate loci, even with equal IDs.
        return self.record_index, self.feature_type, self.location_key


def location_key(location):
    # Preserve biological part order and per-part strand for trans-spliced genes.
    return (
        getattr(location, "operator", None),
        tuple((str(p.start), str(p.end), p.strand, p.ref, p.ref_db)
              for p in location.parts),
    )


def parse_features(path):
    require_biopython()
    features = []
    with StringIO(normalized_genbank_text_for_biopython(path)) as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found: {path}")
    for record_index, record in enumerate(records):
        for feature_index, feature in enumerate(record.features):
            if feature.type not in DUPLICATE_TYPES:
                continue
            if feature.location is None:
                raise ValueError(
                    f"Unparseable location in {path}, record {record.id}, "
                    f"feature {feature_index} ({feature.type})"
                )
            q = feature.qualifiers
            gene = str((q.get("gene") or [""])[0]).strip()
            rendered = feature_location_string(feature, len(record.seq))
            features.append(Feature(
                record_id=record_locus_id(record), record_index=record_index,
                feature_index=feature_index, gene=gene, feature_type=feature.type,
                location=rendered, location_key=location_key(feature.location),
                intervals=tuple((int(p.start), int(p.end)) for p in feature.location.parts),
                locus_tag=str((q.get("locus_tag") or [""])[0]),
                pseudogene="pseudo" in q or "pseudogene" in q,
                partial="<" in rendered or ">" in rendered,
            ))
    return features


def inventory_features(features):
    products = [f for f in features if f.feature_type in PRODUCT_TYPES and f.gene]
    # A gene-only pseudogene is still annotated, rather than a missing gene.
    for gene in features:
        if gene.feature_type != "gene" or not gene.pseudogene or not gene.gene:
            continue
        represented = any(
            f.feature_type in PRODUCT_TYPES
            and f.record_index == gene.record_index and f.gene_key == gene.gene_key
            and (
                (gene.locus_tag and gene.locus_tag == f.locus_tag)
                or gene.location_key == f.location_key
                or all(any(p[0] >= g[0] and p[1] <= g[1]
                           for g in gene.intervals) for p in f.intervals)
            )
            for f in products
        )
        if not represented:
            products.append(gene)
    grouped = defaultdict(list)
    for feature in products:
        grouped[feature.gene_key].append(feature)
    return grouped


def unique_copies(features):
    return {f.copy_key for f in features}


def duplicate_features(features):
    groups = defaultdict(list)
    for feature in features:
        groups[(feature.record_index, feature.gene_key,
                feature.feature_type, feature.location_key)].append(feature)
    rows = []
    for group in groups.values():
        if len(group) < 2:
            continue
        f = group[0]
        rows.append({
            "record_id": f.record_id, "gene": f.gene,
            "feature_type": f.feature_type, "location": f.location,
            "feature_indices": [x.feature_index for x in group],
            "locus_tags": [x.locus_tag for x in group],
            "feature_count": len(group),
        })
    return rows


def reference_files(reference_dir):
    if reference_dir is None:
        return []
    if not reference_dir.is_dir():
        raise ValueError(f"Reference directory does not exist: {reference_dir}")
    paths = sorted(p for p in reference_dir.iterdir()
                   if p.is_file() and p.suffix.lower() in REFERENCE_EXTENSIONS)
    if not paths:
        raise ValueError(f"No GenBank reference files found: {reference_dir}")
    return paths


def build_inventory(annotation, references):
    features = parse_features(annotation)
    target = inventory_features(features)
    ref_groups = [(str(path), inventory_features(parse_features(path)))
                  for path in references]
    duplicates = duplicate_features(features)
    duplicate_genes = {r["gene"].casefold() for r in duplicates if r["gene"]}
    gene_keys = set(target)
    for _, group in ref_groups:
        gene_keys.update(group)
    rows = []
    for key in sorted(gene_keys):
        found = target.get(key, [])
        ref_found = [f for _, group in ref_groups for f in group.get(key, [])]
        gene = (found or ref_found)[0].gene
        count = len(unique_copies(found))
        counts = [len(unique_copies(group.get(key, []))) for _, group in ref_groups]
        present = sum(c > 0 for c in counts)
        reasons = []
        if present and not count:
            reasons.append("missing_candidate")
        elif present and not min(counts) <= count <= max(counts):
            reasons.append("copy_number_difference")
        if key in duplicate_genes:
            reasons.append("duplicate_annotation")
        rows.append({
            "gene": gene,
            "target_feature_types": sorted({f.feature_type for f in found}),
            "target_copy_count": count,
            "target_feature_count": len(found),
            "target_pseudogene_count": len(unique_copies([f for f in found if f.pseudogene])),
            "target_partial_count": len(unique_copies([f for f in found if f.partial])),
            "reference_count": len(references),
            "reference_present_count": present,
            "reference_copy_min": min(counts) if counts else None,
            "reference_copy_max": max(counts) if counts else None,
            "reference_copy_counts": dict(zip((p for p, _ in ref_groups), counts)),
            "target_locations": [{
                "record_id": f.record_id, "feature_index": f.feature_index,
                "feature_type": f.feature_type, "location": f.location,
                "locus_tag": f.locus_tag,
            } for f in found],
            "reference_comparison": (
                "not_evaluated" if not references else
                "target_only" if not present else "evaluated"
            ),
            "review_reason": ";".join(reasons),
            "manual_review_candidate": bool(reasons),
        })
    return {
        "schema_version": 1, "annotation": str(annotation),
        "references": [str(p) for p in references],
        "reference_comparison_evaluated": bool(references),
        "unnamed_target_product_features": sum(
            f.feature_type in PRODUCT_TYPES and not f.gene for f in features),
        "missing_candidate_count": sum("missing_candidate" in r["review_reason"] for r in rows),
        "copy_number_difference_count": sum("copy_number_difference" in r["review_reason"] for r in rows),
        "duplicate_group_count": len(duplicates),
        "manual_review_candidate_count": sum(r["manual_review_candidate"] for r in rows),
        "inventory": rows, "duplicates": duplicates,
    }


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, sort_keys=True)
                             if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def run_inventory(*, annotation, reference_dir, inventory_tsv, duplicates_tsv, summary_json):
    result = build_inventory(annotation, reference_files(reference_dir))
    write_tsv(inventory_tsv, INVENTORY_FIELDS, result["inventory"])
    write_tsv(duplicates_tsv, DUPLICATE_FIELDS, result["duplicates"])
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--inventory-tsv", type=Path, required=True)
    parser.add_argument("--duplicates-tsv", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    args = parser.parse_args()
    run_inventory(**vars(args))


if __name__ == "__main__":
    main()
