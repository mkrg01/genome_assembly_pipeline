import argparse
import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from rename_submission_gene_models import make_submission_prefix


QUALIFIER_INDENT = "                     "
TARGET_FEATURE_KEYS = {
    "gene",
    "mRNA",
    "CDS",
    "5'UTR",
    "3'UTR",
    "intron",
    "exon",
    "tRNA",
    "rRNA",
    "ncRNA",
    "tmRNA",
    "precursor_RNA",
    "prim_transcript",
    "misc_RNA",
}
ORGANELLE_CODES = {
    "mitochondrion": "MT",
    "mitochondria": "MT",
    "mito": "MT",
    "chloroplast": "CP",
    "plastid": "CP",
    "pltd": "CP",
}
FEATURE_LINE_RE = re.compile(r"^     (\S+)\s+(.+)")
QUALIFIER_RE_TEMPLATE = r"^\s*/{}=(?:\"([^\"]*)\"|([^ \t\r\n]+))"
RANGE_RE = re.compile(r"<?(\d+)\.\.>?(\d+)")
POSITION_RE = re.compile(r"(?<![A-Za-z_])<?(\d+)(?![A-Za-z_])")


@dataclass
class FeatureBlock:
    key: str
    start: int
    end: int
    lines: list[str]
    record_index: int
    feature_index: int
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=dict)
    intervals: list[tuple[int, int]] = field(default_factory=list)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add DDBJ-style organelle /locus_tag qualifiers to GenBank annotation "
            "features and write the result as plain text or gzip."
        )
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--species-name",
        required=True,
        help=(
            "Scientific name used to derive the registered prefix, "
            'e.g. "Triantha_japonica".'
        ),
    )
    parser.add_argument("--organelle", required=True)
    return parser.parse_args()


def open_text_read(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def open_text_write(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt")
    return path.open("w")


def make_organelle_locus_tag_prefix(species_name: str, organelle: str):
    try:
        organelle_code = ORGANELLE_CODES[organelle]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported organelle '{organelle}'. Expected one of: "
            f"{', '.join(sorted(ORGANELLE_CODES))}"
        ) from exc
    return f"{make_submission_prefix(species_name)}_{organelle_code}"


def qualifier_values(block_lines: list[str], qualifier: str):
    pattern = re.compile(QUALIFIER_RE_TEMPLATE.format(re.escape(qualifier)))
    values = []
    for line in block_lines[1:]:
        stripped = line.strip()
        if not stripped.startswith(f"/{qualifier}="):
            continue
        match = pattern.match(stripped)
        if match:
            values.append(match.group(1) if match.group(1) is not None else match.group(2))
    return values


def feature_location(block_lines: list[str]):
    match = FEATURE_LINE_RE.match(block_lines[0])
    if not match:
        return ""
    parts = [match.group(2).strip()]
    for line in block_lines[1:]:
        stripped = line.strip()
        if stripped.startswith("/"):
            break
        if stripped:
            parts.append(stripped)
    return "".join(parts)


def location_intervals(location: str):
    intervals = [
        (min(int(start), int(end)), max(int(start), int(end)))
        for start, end in RANGE_RE.findall(location)
    ]
    if intervals:
        return intervals

    positions = [int(position) for position in POSITION_RE.findall(location)]
    return [(position, position) for position in positions]


def location_start(intervals: list[tuple[int, int]]):
    if not intervals:
        return sys.maxsize
    return min(start for start, _ in intervals)


def location_end(intervals: list[tuple[int, int]]):
    if not intervals:
        return sys.maxsize
    return max(end for _, end in intervals)


def interval_overlap(left: list[tuple[int, int]], right: list[tuple[int, int]]):
    overlap = 0
    for left_start, left_end in left:
        for right_start, right_end in right:
            start = max(left_start, right_start)
            end = min(left_end, right_end)
            if start <= end:
                overlap += end - start + 1
    return overlap


def parse_feature_blocks(record_lines: list[str], record_index: int):
    in_features = False
    feature_start = None
    feature_index = 0
    blocks = []

    for line_index, line in enumerate(record_lines):
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN") or line.startswith("BASE COUNT") or line.startswith("//"):
            if feature_start is not None:
                blocks.append((feature_start, line_index))
            break
        match = FEATURE_LINE_RE.match(line)
        if not match:
            continue
        if feature_start is not None:
            blocks.append((feature_start, line_index))
        feature_start = line_index

    feature_blocks = []
    for start, end in blocks:
        lines = record_lines[start:end]
        match = FEATURE_LINE_RE.match(lines[0])
        if not match:
            continue
        key = match.group(1)
        location = feature_location(lines)
        feature_blocks.append(
            FeatureBlock(
                key=key,
                start=start,
                end=end,
                lines=lines,
                record_index=record_index,
                feature_index=feature_index,
                location=location,
                qualifiers={
                    "gene": qualifier_values(lines, "gene"),
                    "locus_tag": qualifier_values(lines, "locus_tag"),
                },
                intervals=location_intervals(location),
            )
        )
        feature_index += 1

    return feature_blocks


def split_records(lines: list[str]):
    records = []
    start = 0
    for index, line in enumerate(lines):
        if line.startswith("//"):
            records.append(lines[start : index + 1])
            start = index + 1
    if start < len(lines):
        records.append(lines[start:])
    return records


def match_gene_feature(target: FeatureBlock, gene_features: list[FeatureBlock]):
    if target.key == "gene":
        return target

    target_genes = set(target.qualifiers.get("gene", []))
    candidates = []
    gene_name_matches = []
    for gene_feature in gene_features:
        gene_names = set(gene_feature.qualifiers.get("gene", []))
        if target_genes and gene_names and target_genes.isdisjoint(gene_names):
            continue
        if target_genes and gene_names:
            gene_name_matches.append(gene_feature)
        overlap = interval_overlap(target.intervals, gene_feature.intervals)
        if overlap <= 0:
            continue
        candidates.append((overlap, -gene_feature.feature_index, gene_feature))

    if candidates:
        return max(candidates, key=lambda candidate: (candidate[0], candidate[1]))[2]
    if len(gene_name_matches) == 1:
        return gene_name_matches[0]
    return None


def assign_locus_tags(records: list[list[str]], locus_tag_prefix: str):
    records_by_index = []
    tag_groups = {}
    group_features = {}

    for record_index, record_lines in enumerate(records):
        features = parse_feature_blocks(record_lines, record_index)
        gene_features = [feature for feature in features if feature.key == "gene"]
        targets = [feature for feature in features if feature.key in TARGET_FEATURE_KEYS]
        records_by_index.append((record_lines, features, targets))

        for target in targets:
            gene_feature = match_gene_feature(target, gene_features)
            if gene_feature is None:
                group_key = ("feature", target.record_index, target.feature_index)
                group_feature = target
            else:
                group_key = ("gene", gene_feature.record_index, gene_feature.feature_index)
                group_feature = gene_feature
            tag_groups[(target.record_index, target.feature_index)] = group_key
            group_features.setdefault(group_key, group_feature)

    sorted_group_keys = sorted(
        group_features,
        key=lambda group_key: (
            group_features[group_key].record_index,
            location_start(group_features[group_key].intervals),
            location_end(group_features[group_key].intervals),
            group_features[group_key].feature_index,
        ),
    )
    locus_tags = {
        group_key: f"{locus_tag_prefix}{index:06d}"
        for index, group_key in enumerate(sorted_group_keys, start=1)
    }

    updated_records = []
    changed_feature_count = 0
    for record_lines, features, targets in records_by_index:
        replacement_by_start = {}
        for target in targets:
            target_key = (target.record_index, target.feature_index)
            locus_tag = locus_tags[tag_groups[target_key]]
            replacement_by_start[target.start] = add_or_replace_locus_tag(
                target.lines,
                locus_tag,
            )
            changed_feature_count += 1

        updated_lines = []
        cursor = 0
        for feature in features:
            updated_lines.extend(record_lines[cursor : feature.start])
            if feature.start in replacement_by_start:
                updated_lines.extend(replacement_by_start[feature.start])
            else:
                updated_lines.extend(record_lines[feature.start : feature.end])
            cursor = feature.end
        updated_lines.extend(record_lines[cursor:])
        updated_records.append(updated_lines)

    return [line for record in updated_records for line in record], {
        "locus_tag_prefix": locus_tag_prefix,
        "gene_count": len(sorted_group_keys),
        "feature_count": changed_feature_count,
    }


def add_or_replace_locus_tag(block_lines: list[str], locus_tag: str):
    locus_tag_line = f'{QUALIFIER_INDENT}/locus_tag="{locus_tag}"\n'
    without_locus_tag = [
        line for line in block_lines if not line.strip().startswith("/locus_tag=")
    ]
    insert_index = None
    for index, line in enumerate(without_locus_tag):
        if line.strip().startswith("/gene="):
            insert_index = index + 1
            break
    if insert_index is None:
        insert_index = 1
        while insert_index < len(without_locus_tag):
            stripped = without_locus_tag[insert_index].strip()
            if stripped.startswith("/"):
                break
            insert_index += 1
    return (
        without_locus_tag[:insert_index]
        + [locus_tag_line]
        + without_locus_tag[insert_index:]
    )


def add_organelle_locus_tags(
    input_path: Path,
    output_path: Path,
    species_name: str,
    organelle: str,
):
    locus_tag_prefix = make_organelle_locus_tag_prefix(species_name, organelle)
    with open_text_read(input_path) as handle:
        lines = handle.readlines()

    updated_lines, summary = assign_locus_tags(split_records(lines), locus_tag_prefix)

    with open_text_write(output_path) as handle:
        handle.writelines(updated_lines)

    return summary


def main():
    args = parse_args()
    summary = add_organelle_locus_tags(
        args.input,
        args.output,
        args.species_name,
        args.organelle,
    )
    print(
        f"Added {summary['locus_tag_prefix']} locus_tags to "
        f"{summary['feature_count']} features across {summary['gene_count']} genes.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
