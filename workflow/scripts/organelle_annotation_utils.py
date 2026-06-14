import json
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path


GENBANK_SUFFIXES = (".gb", ".gbk", ".gbf")
FASTA_LINE_WIDTH = 80
FEATURE_INDENT = "     "
QUALIFIER_INDENT = "                     "
SPECIES_PLACEHOLDER = "YOUR_SPECIES"
GENBANK_TOPOLOGIES = {"circular", "linear"}
FASTA_CIRCULAR_PATTERN = re.compile(r"(?:^|\s)circular=(true|false)(?=\s|$)", re.I)
GENBANK_DATE_PATTERN = re.compile(r"\d{2}-[A-Z]{3}-\d{4}")
GENBANK_MOLECULE_PATTERN = re.compile(r"\b(mRNA|DNA|RNA)\b")
GENBANK_FEATURE_LINE_PATTERN = re.compile(r"^     (\S+)\s+(.+)")
GENBANK_QUALIFIER_PATTERN = re.compile(r"^/([^=\s]+)(?:=|$)")
GENBANK_LOCATION_RANGE_PATTERN = re.compile(r"<?(\d+)\.\.>?(\d+)")
GENBANK_LOCATION_RANGE_DETAIL_PATTERN = re.compile(r"(<)?(\d+)\.\.(>)?(\d+)")
GENBANK_LOCATION_POSITION_PATTERN = re.compile(r"(?<![A-Za-z_])<?(\d+)(?![A-Za-z_])")
GENBANK_FEATURE_LOCATION_WIDTH = 58
FEATURE_SORT_ORDER = {
    "gene": 0,
    "mRNA": 1,
    "CDS": 1,
    "tRNA": 1,
    "rRNA": 1,
    "ncRNA": 1,
    "tmRNA": 1,
    "precursor_RNA": 1,
    "prim_transcript": 1,
    "misc_RNA": 1,
    "exon": 2,
    "5'UTR": 2,
    "3'UTR": 2,
    "intron": 2,
}
STANDARD_GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
GENETIC_CODES = {
    "1": STANDARD_GENETIC_CODE,
    "11": STANDARD_GENETIC_CODE,
}
GENETIC_CODE_START_CODONS = {
    "1": {"TTG", "CTG", "ATG"},
    "11": {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
}
GENETIC_CODE_STOP_CODONS = {
    "1": {"TAA", "TAG", "TGA"},
    "11": {"TAA", "TAG", "TGA"},
}
DNA_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
)


@dataclass
class GenBankFeatureBlock:
    key: str
    start: int
    end: int
    lines: list[str]
    feature_index: int
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=dict)
    intervals: list[tuple[int, int]] = field(default_factory=list)


def fasta_record_id(header: str):
    record_id = header.split()[0] if header.split() else ""
    if not record_id:
        raise ValueError(f"Could not determine a FASTA record ID from: {header!r}")
    return record_id


def make_fasta_record(header: str, sequence_lines: list[str]):
    return {
        "header": header,
        "id": fasta_record_id(header),
        "sequence": "".join(sequence_lines).upper().replace("U", "T"),
        "topology": topology_from_fasta_header(header),
    }


def parse_fasta_records(path: Path):
    records = []
    header = None
    sequence_lines = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(make_fasta_record(header, sequence_lines))
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Empty FASTA header found in {path} at line {line_number}."
                    )
                sequence_lines = []
                continue

            if header is None:
                if line.strip():
                    raise ValueError(
                        f"Sequence data found before the first FASTA header in {path} "
                        f"at line {line_number}."
                    )
                continue
            sequence_lines.append("".join(line.split()))

    if header is not None:
        records.append(make_fasta_record(header, sequence_lines))
    if not records:
        raise ValueError(f"No FASTA records found in {path}.")
    return records


def write_fasta_records(records: list[dict], path: Path, line_width=FASTA_LINE_WIDTH):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(f">{record['header']}\n")
            sequence = record["sequence"]
            for index in range(0, len(sequence), line_width):
                handle.write(f"{sequence[index:index + line_width]}\n")


def rotate_sequence(sequence: str, offset: int):
    if not sequence:
        return sequence
    offset %= len(sequence)
    if offset == 0:
        return sequence
    return sequence[offset:] + sequence[:offset]


def rotate_fasta_record(record: dict, offset: int):
    rotated = dict(record)
    rotated["sequence"] = rotate_sequence(record["sequence"], offset)
    return rotated


def organism_name_from_assembly_name(assembly_name: str):
    organism_name = " ".join(assembly_name.replace("_", " ").split())
    if not organism_name:
        raise ValueError("Assembly name must not be empty when deriving organism name.")
    if any(character in organism_name for character in '"\r\n'):
        raise ValueError(
            f"Assembly name cannot be used as a GenBank organism: {assembly_name!r}"
        )
    return organism_name


def copy_first_genbank(output_dir: Path, annotation_path: Path):
    candidates = sorted(
        path
        for path in output_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES
    )
    if not candidates:
        raise RuntimeError(f"No GenBank output was found under {output_dir}.")

    annotation_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(candidates[0], annotation_path)
    return candidates[0]


def split_genbank_records(lines):
    records = []
    current = []
    for line in lines:
        current.append(line)
        if line.strip() == "//":
            records.append(current)
            current = []
    if current:
        records.append(current)
    return records


def trim_genbank_record_to_core_sections(record, path: Path):
    locus_index = None
    features_index = None
    origin_index = None
    for index, line in enumerate(record):
        if line.startswith("LOCUS") and locus_index is None:
            locus_index = index
        elif line.startswith("FEATURES") and features_index is None:
            features_index = index
        elif line.startswith("ORIGIN") and origin_index is None:
            origin_index = index

    missing = [
        name
        for name, index in (
            ("LOCUS", locus_index),
            ("FEATURES", features_index),
            ("ORIGIN", origin_index),
        )
        if index is None
    ]
    if missing:
        raise ValueError(
            f"GenBank record in {path} is missing required section(s): "
            + ", ".join(missing)
        )
    if not locus_index < features_index < origin_index:
        raise ValueError(
            f"GenBank record in {path} has sections out of order: "
            "expected LOCUS, FEATURES, ORIGIN."
        )

    return [record[locus_index], *record[features_index:]]


def trim_genbank_to_core_sections(path: Path):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    trimmed_records = [
        trim_genbank_record_to_core_sections(record, path) for record in records
    ]
    trimmed_lines = [line for record in trimmed_records for line in record]
    changed = lines != trimmed_lines
    if changed:
        path.write_text("".join(trimmed_lines))
    return {
        "core_sections": ["LOCUS", "FEATURES", "ORIGIN"],
        "core_sections_record_count": len(trimmed_records),
        "core_sections_changed": changed,
    }


def genbank_qualifier_line_content(line: str):
    raw = line.rstrip("\r\n")
    if raw.startswith(QUALIFIER_INDENT):
        return raw[len(QUALIFIER_INDENT) :]
    return raw.strip()


def genbank_qualifier_blocks(block_lines: list[str]):
    blocks = []
    current_name = None
    current_lines = []
    for line in block_lines[1:]:
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is not None:
            if current_name is not None:
                blocks.append((current_name, current_lines))
            current_name = qualifier_name
            current_lines = [line]
            continue
        if current_name is not None:
            current_lines.append(line)
    if current_name is not None:
        blocks.append((current_name, current_lines))
    return blocks


def genbank_qualifier_block_value(qualifier: str, lines: list[str]):
    if not lines:
        return None
    first_content = genbank_qualifier_line_content(lines[0])
    prefix = f"/{qualifier}"
    if first_content == prefix:
        return None
    if not first_content.startswith(f"{prefix}="):
        return None
    value_parts = [first_content.removeprefix(f"{prefix}=")]
    value_parts.extend(genbank_qualifier_line_content(line) for line in lines[1:])
    value = "".join(value_parts)
    if value.startswith('"'):
        value = value[1:]
        if value.endswith('"'):
            value = value[:-1]
    return value


def genbank_qualifier_values(block_lines: list[str], qualifier: str):
    values = []
    for qualifier_name, qualifier_lines in genbank_qualifier_blocks(block_lines):
        if qualifier_name != qualifier:
            continue
        value = genbank_qualifier_block_value(qualifier, qualifier_lines)
        if value is not None:
            values.append(value)
    return values


def genbank_feature_location(block_lines: list[str]):
    match = GENBANK_FEATURE_LINE_PATTERN.match(block_lines[0])
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


def genbank_location_intervals(location: str):
    intervals = [
        (min(int(start), int(end)), max(int(start), int(end)))
        for start, end in GENBANK_LOCATION_RANGE_PATTERN.findall(location)
    ]
    if intervals:
        return intervals

    positions = [
        int(position) for position in GENBANK_LOCATION_POSITION_PATTERN.findall(location)
    ]
    return [(position, position) for position in positions]


def genbank_location_start(intervals: list[tuple[int, int]]):
    if not intervals:
        return float("inf")
    return min(start for start, _ in intervals)


def genbank_location_end(intervals: list[tuple[int, int]]):
    if not intervals:
        return float("inf")
    return max(end for _, end in intervals)


def genbank_interval_overlap(left: list[tuple[int, int]], right: list[tuple[int, int]]):
    overlap = 0
    for left_start, left_end in left:
        for right_start, right_end in right:
            start = max(left_start, right_start)
            end = min(left_end, right_end)
            if start <= end:
                overlap += end - start + 1
    return overlap


def parse_genbank_feature_blocks(record: list[str]):
    features_index = None
    origin_index = None
    feature_start = None
    block_ranges = []

    for index, line in enumerate(record):
        if line.startswith("FEATURES") and features_index is None:
            features_index = index
            continue
        if features_index is None:
            continue
        if line.startswith("ORIGIN") or line.startswith("BASE COUNT") or line.startswith("//"):
            origin_index = index
            if feature_start is not None:
                block_ranges.append((feature_start, index))
            break
        if GENBANK_FEATURE_LINE_PATTERN.match(line):
            if feature_start is not None:
                block_ranges.append((feature_start, index))
            feature_start = index

    if features_index is None or origin_index is None:
        return features_index, origin_index, []

    blocks = []
    for feature_index, (start, end) in enumerate(block_ranges):
        lines = record[start:end]
        match = GENBANK_FEATURE_LINE_PATTERN.match(lines[0])
        if match is None:
            continue
        location = genbank_feature_location(lines)
        blocks.append(
            GenBankFeatureBlock(
                key=match.group(1),
                start=start,
                end=end,
                lines=lines,
                feature_index=feature_index,
                location=location,
                qualifiers={"gene": genbank_qualifier_values(lines, "gene")},
                intervals=genbank_location_intervals(location),
            )
        )

    return features_index, origin_index, blocks


def genbank_feature_qualifier_start(block_lines: list[str]):
    for index, line in enumerate(block_lines[1:], start=1):
        if line.strip().startswith("/"):
            return index
    return len(block_lines)


def split_genbank_location_at_commas(location: str):
    parts = []
    start = 0
    for index, character in enumerate(location):
        if character != ",":
            continue
        parts.append(location[start : index + 1])
        start = index + 1
    if start < len(location):
        parts.append(location[start:])
    return [part for part in parts if part]


def wrap_genbank_feature_location(
    location: str,
    width: int = GENBANK_FEATURE_LOCATION_WIDTH,
):
    parts = split_genbank_location_at_commas(location)
    if not parts:
        return [location]

    wrapped = []
    current = ""
    for part in parts:
        if not current:
            current = part
            continue
        if len(current) + len(part) <= width:
            current += part
            continue
        wrapped.append(current)
        current = part

    if current:
        wrapped.append(current)
    return wrapped


def format_genbank_feature_location_lines(
    key: str,
    location: str,
    newline: str = "\n",
):
    location_parts = wrap_genbank_feature_location(location)
    first_prefix = f"{FEATURE_INDENT}{key:<16}"
    return [
        f"{first_prefix if index == 0 else QUALIFIER_INDENT}{part}{newline}"
        for index, part in enumerate(location_parts)
    ]


def normalize_genbank_record_feature_location_wrapping(record: list[str]):
    features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
    if features_index is None or origin_index is None or not blocks:
        return record, False, 0

    changed_feature_count = 0
    normalized_blocks = []
    for block in blocks:
        qualifier_start = genbank_feature_qualifier_start(block.lines)
        newline = "\n" if block.lines[0].endswith("\n") else ""
        normalized_location_lines = format_genbank_feature_location_lines(
            block.key,
            block.location,
            newline=newline,
        )
        if normalized_location_lines != block.lines[:qualifier_start]:
            changed_feature_count += 1
        normalized_blocks.append(
            [*normalized_location_lines, *block.lines[qualifier_start:]]
        )

    if changed_feature_count == 0:
        return record, False, 0

    normalized_feature_lines = [
        line for block_lines in normalized_blocks for line in block_lines
    ]
    normalized_record = [
        *record[: features_index + 1],
        *normalized_feature_lines,
        *record[origin_index:],
    ]
    return normalized_record, True, changed_feature_count


def normalize_genbank_feature_location_wrapping(path: Path):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    normalized_records = []
    changed_record_count = 0
    changed_feature_count = 0

    for record in records:
        normalized_record, changed, record_changed_feature_count = (
            normalize_genbank_record_feature_location_wrapping(record)
        )
        normalized_records.append(normalized_record)
        changed_record_count += int(changed)
        changed_feature_count += record_changed_feature_count

    normalized_lines = [line for record in normalized_records for line in record]
    changed = lines != normalized_lines
    if changed:
        path.write_text("".join(normalized_lines))

    return {
        "feature_location_wrap_record_count": len(records),
        "feature_location_wrap_changed_record_count": changed_record_count,
        "feature_location_wrap_changed_feature_count": changed_feature_count,
        "feature_location_wrap_changed": changed,
    }


def genbank_qualifier_name(line: str):
    match = GENBANK_QUALIFIER_PATTERN.match(line.strip())
    if match is None:
        return None
    return match.group(1)


def remove_genbank_feature_block_qualifiers(
    block_lines: list[str],
    qualifiers: tuple[str, ...],
):
    qualifier_set = set(qualifiers)
    removed_counts = {qualifier: 0 for qualifier in qualifiers}
    cleaned_lines = []
    skipping_qualifier = False

    for line in block_lines:
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is not None:
            if qualifier_name in qualifier_set:
                removed_counts[qualifier_name] += 1
                skipping_qualifier = True
                continue
            skipping_qualifier = False
            cleaned_lines.append(line)
            continue

        if skipping_qualifier:
            continue
        cleaned_lines.append(line)

    return cleaned_lines, removed_counts


def remove_genbank_record_feature_qualifiers(
    record: list[str],
    *,
    feature_keys: tuple[str, ...],
    qualifiers: tuple[str, ...],
):
    features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
    if features_index is None or origin_index is None or not blocks:
        return record, False, 0, 0, {qualifier: 0 for qualifier in qualifiers}

    feature_key_set = set(feature_keys)
    target_feature_count = 0
    changed_feature_count = 0
    removed_counts = {qualifier: 0 for qualifier in qualifiers}
    cleaned_blocks = []

    for block in blocks:
        if block.key not in feature_key_set:
            cleaned_blocks.append(block.lines)
            continue

        target_feature_count += 1
        cleaned_lines, block_removed_counts = remove_genbank_feature_block_qualifiers(
            block.lines,
            qualifiers,
        )
        if cleaned_lines != block.lines:
            changed_feature_count += 1
        for qualifier, count in block_removed_counts.items():
            removed_counts[qualifier] += count
        cleaned_blocks.append(cleaned_lines)

    if changed_feature_count == 0:
        return record, False, target_feature_count, 0, removed_counts

    cleaned_feature_lines = [line for block in cleaned_blocks for line in block]
    cleaned_record = [
        *record[: features_index + 1],
        *cleaned_feature_lines,
        *record[origin_index:],
    ]
    return (
        cleaned_record,
        True,
        target_feature_count,
        changed_feature_count,
        removed_counts,
    )


def remove_genbank_feature_qualifiers(
    path: Path,
    *,
    feature_keys: tuple[str, ...],
    qualifiers: tuple[str, ...],
):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    cleaned_records = []
    changed_records = 0
    target_feature_count = 0
    changed_feature_count = 0
    removed_counts = {qualifier: 0 for qualifier in qualifiers}

    for record in records:
        (
            cleaned_record,
            changed,
            record_target_feature_count,
            record_changed_feature_count,
            record_removed_counts,
        ) = remove_genbank_record_feature_qualifiers(
            record,
            feature_keys=feature_keys,
            qualifiers=qualifiers,
        )
        cleaned_records.append(cleaned_record)
        changed_records += int(changed)
        target_feature_count += record_target_feature_count
        changed_feature_count += record_changed_feature_count
        for qualifier, count in record_removed_counts.items():
            removed_counts[qualifier] += count

    cleaned_lines = [line for record in cleaned_records for line in record]
    changed = lines != cleaned_lines
    if changed:
        path.write_text("".join(cleaned_lines))

    return {
        "qualifier_removal_record_count": len(records),
        "qualifier_removal_changed_record_count": changed_records,
        "qualifier_removal_feature_keys": list(feature_keys),
        "qualifier_removal_qualifiers": list(qualifiers),
        "qualifier_removal_target_feature_count": target_feature_count,
        "qualifier_removal_changed_feature_count": changed_feature_count,
        "qualifier_removal_counts": removed_counts,
        "qualifier_removal_changed": changed,
    }


def genbank_feature_has_qualifier(block_lines: list[str], qualifier: str):
    return any(genbank_qualifier_name(line) == qualifier for line in block_lines[1:])


def normalize_pmga_trans_splicing_feature_block(block_lines: list[str]):
    cleaned_lines = []
    removed_exception_count = 0
    skipping_qualifier = False

    for line in block_lines:
        qualifier_name = genbank_qualifier_name(line)
        if qualifier_name is not None:
            skipping_qualifier = False
            if qualifier_name == "exception":
                values = genbank_qualifier_values([block_lines[0], line], "exception")
                if values and values[0].lower() == "trans-splicing":
                    removed_exception_count += 1
                    skipping_qualifier = True
                    continue
            cleaned_lines.append(line)
            continue

        if skipping_qualifier:
            continue
        cleaned_lines.append(line)

    added_trans_splicing = False
    if removed_exception_count and not genbank_feature_has_qualifier(
        cleaned_lines,
        "trans_splicing",
    ):
        cleaned_lines.append(f"{QUALIFIER_INDENT}/trans_splicing\n")
        added_trans_splicing = True

    return cleaned_lines, removed_exception_count, added_trans_splicing


def normalize_pmga_trans_splicing_record(
    record: list[str],
    *,
    feature_keys: tuple[str, ...] = ("CDS",),
):
    features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
    if features_index is None or origin_index is None or not blocks:
        return record, False, 0, 0, 0

    feature_key_set = set(feature_keys)
    changed_feature_count = 0
    removed_exception_count = 0
    added_trans_splicing_count = 0
    cleaned_blocks = []

    for block in blocks:
        if block.key not in feature_key_set:
            cleaned_blocks.append(block.lines)
            continue

        cleaned_lines, block_removed_count, block_added = (
            normalize_pmga_trans_splicing_feature_block(block.lines)
        )
        if cleaned_lines != block.lines:
            changed_feature_count += 1
        removed_exception_count += block_removed_count
        added_trans_splicing_count += int(block_added)
        cleaned_blocks.append(cleaned_lines)

    if changed_feature_count == 0:
        return record, False, 0, removed_exception_count, added_trans_splicing_count

    cleaned_feature_lines = [line for block in cleaned_blocks for line in block]
    cleaned_record = [
        *record[: features_index + 1],
        *cleaned_feature_lines,
        *record[origin_index:],
    ]
    return (
        cleaned_record,
        True,
        changed_feature_count,
        removed_exception_count,
        added_trans_splicing_count,
    )


def normalize_pmga_trans_splicing_qualifiers(
    path: Path,
    *,
    feature_keys: tuple[str, ...] = ("CDS",),
):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    cleaned_records = []
    changed_record_count = 0
    changed_feature_count = 0
    removed_exception_count = 0
    added_trans_splicing_count = 0

    for record in records:
        (
            cleaned_record,
            changed,
            record_changed_feature_count,
            record_removed_exception_count,
            record_added_trans_splicing_count,
        ) = normalize_pmga_trans_splicing_record(
            record,
            feature_keys=feature_keys,
        )
        cleaned_records.append(cleaned_record)
        changed_record_count += int(changed)
        changed_feature_count += record_changed_feature_count
        removed_exception_count += record_removed_exception_count
        added_trans_splicing_count += record_added_trans_splicing_count

    cleaned_lines = [line for record in cleaned_records for line in record]
    changed = lines != cleaned_lines
    if changed:
        path.write_text("".join(cleaned_lines))

    return {
        "pmga_trans_splicing_record_count": len(records),
        "pmga_trans_splicing_changed_record_count": changed_record_count,
        "pmga_trans_splicing_feature_keys": list(feature_keys),
        "pmga_trans_splicing_changed_feature_count": changed_feature_count,
        "pmga_trans_splicing_removed_exception_count": removed_exception_count,
        "pmga_trans_splicing_added_qualifier_count": added_trans_splicing_count,
        "pmga_trans_splicing_changed": changed,
    }


def first_genbank_qualifier_value(
    block_lines: list[str],
    qualifier: str,
    default=None,
):
    values = genbank_qualifier_values(block_lines, qualifier)
    if not values:
        return default
    return values[0]


def genbank_record_origin_sequence(record: list[str]):
    origin_index = None
    for index, line in enumerate(record):
        if line.startswith("ORIGIN"):
            origin_index = index
            break
    if origin_index is None:
        return ""

    sequence_parts = []
    for line in record[origin_index + 1 :]:
        if line.strip() == "//":
            break
        sequence_parts.extend(re.findall(r"[A-Za-z]+", line))
    return "".join(sequence_parts).upper().replace("U", "T")


def reverse_complement_dna(sequence: str):
    return sequence.translate(DNA_COMPLEMENT_TABLE)[::-1]


def unwrap_genbank_location_function(location: str):
    location = location.strip()
    for function_name in ("complement", "join", "order"):
        prefix = f"{function_name}("
        if not location.startswith(prefix) or not location.endswith(")"):
            continue
        depth = 0
        for index, character in enumerate(location[len(function_name) :], start=len(function_name)):
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0 and index != len(location) - 1:
                    break
        else:
            if depth == 0:
                return function_name, location[len(prefix) : -1]
    return None, None


def split_genbank_location_arguments(location: str):
    parts = []
    start = 0
    depth = 0
    for index, character in enumerate(location):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif character == "," and depth == 0:
            parts.append(location[start:index].strip())
            start = index + 1
    parts.append(location[start:].strip())
    return [part for part in parts if part]


def clean_genbank_location_position(position: str):
    return int(position.strip().lstrip("<>"))


def extract_genbank_location_sequence(location: str, record_sequence: str):
    normalized = "".join(location.split())
    function_name, inner = unwrap_genbank_location_function(normalized)
    if function_name == "complement":
        return reverse_complement_dna(
            extract_genbank_location_sequence(inner, record_sequence)
        )
    if function_name in {"join", "order"}:
        return "".join(
            extract_genbank_location_sequence(part, record_sequence)
            for part in split_genbank_location_arguments(inner)
        )

    if "," in normalized:
        return "".join(
            extract_genbank_location_sequence(part, record_sequence)
            for part in split_genbank_location_arguments(normalized)
        )

    if ":" in normalized or "^" in normalized:
        raise ValueError(f"Unsupported GenBank location: {location}")

    range_match = re.fullmatch(r"([<>]?\d+)\.\.([<>]?\d+)", normalized)
    if range_match:
        start = clean_genbank_location_position(range_match.group(1))
        end = clean_genbank_location_position(range_match.group(2))
        if start < 1 or end < 1 or start > end:
            raise ValueError(f"Unsupported GenBank location range: {location}")
        return record_sequence[start - 1 : end]

    position_match = re.fullmatch(r"[<>]?(\d+)", normalized)
    if position_match:
        position = int(position_match.group(1))
        if position < 1:
            raise ValueError(f"Unsupported GenBank location position: {location}")
        return record_sequence[position - 1 : position]

    raise ValueError(f"Unsupported GenBank location: {location}")


def translate_complete_codons(sequence: str, table_id: str):
    code = GENETIC_CODES[str(table_id)]
    normalized = sequence.upper().replace("U", "T")
    protein = []
    ambiguous_codons = []
    for index in range(0, len(normalized) - len(normalized) % 3, 3):
        codon = normalized[index : index + 3]
        amino_acid = code.get(codon)
        if amino_acid is None:
            amino_acid = "X"
            ambiguous_codons.append({"codon": codon, "position": index + 1})
        protein.append(amino_acid)
    return "".join(protein), ambiguous_codons


def validate_complete_cds_translation(sequence: str, table_id: str):
    table_id = str(table_id)
    normalized = sequence.upper().replace("U", "T")
    errors = []
    warnings = []

    if not normalized:
        errors.append("empty_cds")
        return {"valid": False, "errors": errors, "warnings": warnings}

    if len(normalized) % 3 != 0:
        errors.append("length_not_multiple_of_three")

    codons = [
        normalized[index : index + 3]
        for index in range(0, len(normalized) - len(normalized) % 3, 3)
    ]
    if not codons:
        errors.append("no_complete_codon")
        return {"valid": False, "errors": errors, "warnings": warnings}

    start_codon = codons[0]
    if start_codon not in GENETIC_CODE_START_CODONS[table_id]:
        errors.append(f"invalid_start_codon:{start_codon}")

    stop_codon = codons[-1]
    if stop_codon not in GENETIC_CODE_STOP_CODONS[table_id]:
        errors.append(f"missing_terminal_stop:{stop_codon}")

    internal_stop_positions = [
        index + 1
        for index, codon in enumerate(codons[:-1])
        if codon in GENETIC_CODE_STOP_CODONS[table_id]
    ]
    if internal_stop_positions:
        positions = ",".join(str(position) for position in internal_stop_positions)
        errors.append(f"internal_stop_codon:{positions}")

    _, ambiguous_codons = translate_complete_codons(normalized, table_id)
    if ambiguous_codons:
        warnings.append(
            "ambiguous_codon:"
            + ",".join(
                f"{item['codon']}@{item['position']}" for item in ambiguous_codons
            )
        )

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
    }


def make_cds_auto_translation_issue(
    *,
    severity: str,
    code: str,
    message: str,
    record_name: str,
    block: GenBankFeatureBlock,
    gene: str,
    details=None,
):
    issue = {
        "severity": severity,
        "code": code,
        "message": message,
        "record": record_name,
        "feature_index": block.feature_index,
        "gene": gene,
        "location": block.location,
    }
    if details:
        issue["details"] = details
    return issue


def validate_genbank_cds_auto_translation(
    path: Path,
    *,
    strict_table: str = "1",
    comparison_table: str = "11",
):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    issues = []
    cds_count = 0
    checked_cds_count = 0
    strict_pass_count = 0
    skipped_exception_count = 0
    skipped_partial_count = 0
    skipped_pseudo_count = 0
    table_compare_match_count = 0
    table_compare_mismatch_count = 0

    for record_index, record in enumerate(records, start=1):
        metadata = genbank_record_locus_metadata(record)
        record_name = metadata["name"] if metadata is not None else f"record_{record_index}"
        record_sequence = genbank_record_origin_sequence(record)
        _, _, blocks = parse_genbank_feature_blocks(record)

        for block in blocks:
            if block.key != "CDS":
                continue
            cds_count += 1
            gene = first_genbank_qualifier_value(block.lines, "gene", "?")
            exceptions = genbank_qualifier_values(block.lines, "exception")
            is_pseudo = genbank_feature_has_qualifier(
                block.lines,
                "pseudo",
            ) or genbank_feature_has_qualifier(block.lines, "pseudogene")
            is_partial = "<" in block.location or ">" in block.location

            if exceptions:
                skipped_exception_count += 1
                continue
            if is_pseudo:
                skipped_pseudo_count += 1
                continue
            if is_partial:
                skipped_partial_count += 1
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="warning",
                        code="partial_cds_skipped",
                        message=(
                            "CDS has a partial location and was not checked as a "
                            "complete auto-translated CDS."
                        ),
                        record_name=record_name,
                        block=block,
                        gene=gene,
                    )
                )
                continue

            checked_cds_count += 1
            try:
                cds_sequence = extract_genbank_location_sequence(
                    block.location,
                    record_sequence,
                )
            except ValueError as error:
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="error",
                        code="location_parse_error",
                        message=str(error),
                        record_name=record_name,
                        block=block,
                        gene=gene,
                    )
                )
                continue

            codon_start = first_genbank_qualifier_value(
                block.lines,
                "codon_start",
                "1",
            )
            try:
                codon_start_value = int(codon_start)
            except (TypeError, ValueError):
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="error",
                        code="invalid_codon_start",
                        message=f"Invalid /codon_start value: {codon_start!r}",
                        record_name=record_name,
                        block=block,
                        gene=gene,
                    )
                )
                continue
            if codon_start_value not in {1, 2, 3}:
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="error",
                        code="invalid_codon_start",
                        message=f"Invalid /codon_start value: {codon_start!r}",
                        record_name=record_name,
                        block=block,
                        gene=gene,
                    )
                )
                continue

            coding_sequence = cds_sequence[codon_start_value - 1 :]
            strict_translation, strict_ambiguous = translate_complete_codons(
                coding_sequence,
                strict_table,
            )
            comparison_translation, comparison_ambiguous = translate_complete_codons(
                coding_sequence,
                comparison_table,
            )
            if strict_translation == comparison_translation:
                table_compare_match_count += 1
            else:
                table_compare_mismatch_count += 1
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="error",
                        code="translation_table_mismatch",
                        message=(
                            f"Conceptual translation differs between table "
                            f"{strict_table} and table {comparison_table}."
                        ),
                        record_name=record_name,
                        block=block,
                        gene=gene,
                        details={
                            f"table_{strict_table}_translation": strict_translation,
                            f"table_{comparison_table}_translation": (
                                comparison_translation
                            ),
                        },
                    )
                )

            strict_validation = validate_complete_cds_translation(
                coding_sequence,
                strict_table,
            )
            comparison_validation = validate_complete_cds_translation(
                coding_sequence,
                comparison_table,
            )
            if strict_validation["valid"]:
                strict_pass_count += 1
            else:
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="error",
                        code="strict_cds_translation_failed",
                        message=(
                            f"CDS is not valid for complete auto-translation "
                            f"with table {strict_table}."
                        ),
                        record_name=record_name,
                        block=block,
                        gene=gene,
                        details={
                            f"table_{strict_table}_errors": (
                                strict_validation["errors"]
                            ),
                            f"table_{comparison_table}_valid": (
                                comparison_validation["valid"]
                            ),
                            f"table_{comparison_table}_errors": (
                                comparison_validation["errors"]
                            ),
                        },
                    )
                )

            ambiguous = strict_ambiguous or comparison_ambiguous
            if ambiguous or strict_validation["warnings"]:
                issues.append(
                    make_cds_auto_translation_issue(
                        severity="warning",
                        code="ambiguous_codon",
                        message="CDS contains ambiguous codon(s) translated as X.",
                        record_name=record_name,
                        block=block,
                        gene=gene,
                        details={
                            "ambiguous_codons": strict_ambiguous,
                            "warnings": strict_validation["warnings"],
                        },
                    )
                )

    blocking_issues = [issue for issue in issues if issue["severity"] == "error"]
    warning_issues = [issue for issue in issues if issue["severity"] == "warning"]
    return {
        "cds_auto_translation_record_count": len(records),
        "cds_auto_translation_cds_count": cds_count,
        "cds_auto_translation_checked_cds_count": checked_cds_count,
        "cds_auto_translation_strict_table": str(strict_table),
        "cds_auto_translation_comparison_table": str(comparison_table),
        "cds_auto_translation_strict_pass_count": strict_pass_count,
        "cds_auto_translation_table_compare_match_count": table_compare_match_count,
        "cds_auto_translation_table_compare_mismatch_count": (
            table_compare_mismatch_count
        ),
        "cds_auto_translation_skipped_exception_count": skipped_exception_count,
        "cds_auto_translation_skipped_partial_count": skipped_partial_count,
        "cds_auto_translation_skipped_pseudo_count": skipped_pseudo_count,
        "cds_auto_translation_blocking_issue_count": len(blocking_issues),
        "cds_auto_translation_warning_count": len(warning_issues),
        "cds_auto_translation_issues": issues,
        "cds_auto_translation_passed": not blocking_issues,
    }


def format_cds_auto_translation_validation_errors(validation):
    blocking_issues = [
        issue
        for issue in validation.get("cds_auto_translation_issues", [])
        if issue.get("severity") == "error"
    ]
    issue_lines = []
    for issue in blocking_issues[:10]:
        issue_lines.append(
            "- "
            f"{issue.get('record')} {issue.get('gene')} "
            f"({issue.get('code')}): {issue.get('message')} "
            f"location={issue.get('location')}"
        )
    if len(blocking_issues) > 10:
        issue_lines.append(f"- ... {len(blocking_issues) - 10} more issue(s)")
    joined_issues = "\n".join(issue_lines)
    return (
        "CDS auto-translation validation failed before removing "
        "/translation and /transl_table from PMGA CDS features.\n"
        f"{joined_issues}"
    )


def normalize_circular_origin_wrapping_location(
    location: str,
    sequence_length: int,
):
    replacements = 0

    def replace_range(match):
        nonlocal replacements
        start_fuzzy, start_text, end_fuzzy, end_text = match.groups()
        start = int(start_text)
        end = int(end_text)
        if start <= end or start > sequence_length or end < 1:
            return match.group(0)

        replacements += 1
        start_prefix = start_fuzzy or ""
        end_prefix = end_fuzzy or ""
        return (
            f"join({start_prefix}{start}..{sequence_length},"
            f"1..{end_prefix}{end})"
        )

    normalized = GENBANK_LOCATION_RANGE_DETAIL_PATTERN.sub(replace_range, location)
    return normalized, replacements


def genbank_record_locus_metadata(record: list[str]):
    for line in record:
        if line.startswith("LOCUS"):
            return extract_locus_metadata(line.rstrip("\n"))
    return None


def clamp_genbank_intervals(
    intervals: list[tuple[int, int]],
    sequence_length: int,
):
    clamped = []
    for start, end in intervals:
        start = max(1, min(sequence_length, start))
        end = max(1, min(sequence_length, end))
        if start <= end:
            clamped.append((start, end))
    return clamped


def merge_genbank_intervals(intervals: list[tuple[int, int]]):
    if not intervals:
        return []

    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
            continue
        merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def genbank_record_non_source_feature_intervals(
    record: list[str],
    sequence_length: int,
):
    _features_index, _origin_index, blocks = parse_genbank_feature_blocks(record)
    intervals = []
    feature_count = 0
    for block in blocks:
        if block.key == "source":
            continue
        clamped = clamp_genbank_intervals(block.intervals, sequence_length)
        if not clamped:
            continue
        feature_count += 1
        intervals.extend(clamped)
    return intervals, feature_count


def select_circular_origin_rotation_from_record(record: list[str]):
    metadata = genbank_record_locus_metadata(record)
    if metadata is None:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "skipped_reason": "no LOCUS line found in preliminary annotation",
        }

    sequence_length = metadata.get("length")
    if metadata.get("topology") != "circular":
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "sequence_length": sequence_length,
            "skipped_reason": "record is not circular",
        }
    if not sequence_length:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "skipped_reason": "could not determine circular sequence length",
        }

    intervals, feature_count = genbank_record_non_source_feature_intervals(
        record,
        sequence_length,
    )
    merged_intervals = merge_genbank_intervals(intervals)
    if len(merged_intervals) < 2:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "sequence_length": sequence_length,
            "feature_count": feature_count,
            "feature_interval_count": len(intervals),
            "merged_feature_interval_count": len(merged_intervals),
            "skipped_reason": (
                "fewer than two non-source feature intervals were available "
                "for internal gap selection"
            ),
        }

    candidates = []
    for (_left_start, left_end), (right_start, _right_end) in zip(
        merged_intervals,
        merged_intervals[1:],
    ):
        gap_start = left_end + 1
        gap_end = right_start - 1
        if gap_start <= gap_end:
            candidates.append(
                {
                    "gap_start": gap_start,
                    "gap_end": gap_end,
                    "gap_length": gap_end - gap_start + 1,
                }
            )

    if not candidates:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "sequence_length": sequence_length,
            "feature_count": feature_count,
            "feature_interval_count": len(intervals),
            "merged_feature_interval_count": len(merged_intervals),
            "candidate_gap_count": 0,
            "skipped_reason": "no internal feature-free gap was found",
        }

    selected = max(candidates, key=lambda item: item["gap_length"])
    new_origin_position = (
        selected["gap_start"] + selected["gap_end"]
    ) // 2
    rotation_offset = new_origin_position - 1
    return {
        "rotation_applied": rotation_offset > 0,
        "rotation_offset": rotation_offset,
        "new_origin_position": new_origin_position,
        "sequence_length": sequence_length,
        "feature_count": feature_count,
        "feature_interval_count": len(intervals),
        "merged_feature_interval_count": len(merged_intervals),
        "candidate_gap_count": len(candidates),
        "selected_gap_start": selected["gap_start"],
        "selected_gap_end": selected["gap_end"],
        "selected_gap_length": selected["gap_length"],
        "skipped_reason": None,
    }


def select_circular_origin_rotation_from_genbank(path: Path):
    records = split_genbank_records(path.read_text().splitlines(keepends=True))
    if not records:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "skipped_reason": "no GenBank records found in preliminary annotation",
        }
    if len(records) > 1:
        return {
            "rotation_applied": False,
            "rotation_offset": 0,
            "new_origin_position": 1,
            "record_count": len(records),
            "skipped_reason": (
                "preliminary annotation contained multiple records; "
                "run circular-origin selection per record"
            ),
        }
    return select_circular_origin_rotation_from_record(records[0])


def normalize_genbank_record_origin_wrapping_locations(record: list[str]):
    features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
    if features_index is None or origin_index is None or not blocks:
        return record, False, 0

    metadata = genbank_record_locus_metadata(record)
    if (
        metadata is None
        or metadata["topology"] != "circular"
        or metadata["length"] is None
    ):
        return record, False, 0

    replacement_count = 0
    normalized_blocks = []
    for block in blocks:
        normalized_location, block_replacements = (
            normalize_circular_origin_wrapping_location(
                block.location,
                metadata["length"],
            )
        )
        replacement_count += block_replacements
        if block_replacements == 0:
            normalized_blocks.append(block.lines)
            continue

        qualifier_start = genbank_feature_qualifier_start(block.lines)
        newline = "\n" if block.lines[0].endswith("\n") else ""
        normalized_blocks.append(
            [
                *format_genbank_feature_location_lines(
                    block.key,
                    normalized_location,
                    newline=newline,
                ),
                *block.lines[qualifier_start:],
            ]
        )

    if replacement_count == 0:
        return record, False, 0

    normalized_feature_lines = [
        line for block_lines in normalized_blocks for line in block_lines
    ]
    normalized_record = [
        *record[: features_index + 1],
        *normalized_feature_lines,
        *record[origin_index:],
    ]
    return normalized_record, True, replacement_count


def normalize_genbank_origin_wrapping_locations(path: Path):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    normalized_records = []
    changed_records = 0
    location_count = 0

    for record in records:
        normalized_record, changed, record_location_count = (
            normalize_genbank_record_origin_wrapping_locations(record)
        )
        normalized_records.append(normalized_record)
        changed_records += int(changed)
        location_count += record_location_count

    normalized_lines = [line for record in normalized_records for line in record]
    changed = lines != normalized_lines
    if changed:
        path.write_text("".join(normalized_lines))

    return {
        "origin_wrap_record_count": len(records),
        "origin_wrap_changed_record_count": changed_records,
        "origin_wrap_location_count": location_count,
        "origin_wrap_changed": changed,
    }


def match_genbank_gene_feature(
    target: GenBankFeatureBlock,
    gene_features: list[GenBankFeatureBlock],
):
    if target.key == "gene":
        return target

    target_genes = set(target.qualifiers.get("gene", []))
    if not target_genes:
        return None

    candidates = []
    gene_name_matches = []
    for gene_feature in gene_features:
        gene_names = set(gene_feature.qualifiers.get("gene", []))
        if target_genes and gene_names and target_genes.isdisjoint(gene_names):
            continue
        if target_genes and gene_names:
            gene_name_matches.append(gene_feature)
        overlap = genbank_interval_overlap(target.intervals, gene_feature.intervals)
        if overlap <= 0:
            continue
        candidates.append((overlap, -gene_feature.feature_index, gene_feature))

    if candidates:
        return max(candidates, key=lambda candidate: (candidate[0], candidate[1]))[2]
    if len(gene_name_matches) == 1:
        return gene_name_matches[0]
    return None


def genbank_feature_sort_key(feature: GenBankFeatureBlock):
    return (
        genbank_location_start(feature.intervals),
        genbank_location_end(feature.intervals),
        feature.feature_index,
    )


def sort_genbank_record_features_by_location(record: list[str]):
    features_index, origin_index, blocks = parse_genbank_feature_blocks(record)
    if features_index is None or origin_index is None or not blocks:
        return record, False, len(blocks), 0

    source_blocks = [block for block in blocks if block.key == "source"]
    sortable_blocks = [block for block in blocks if block.key != "source"]
    gene_features = [block for block in sortable_blocks if block.key == "gene"]
    group_features = {}
    groups = {}

    for block in sortable_blocks:
        gene_feature = match_genbank_gene_feature(block, gene_features)
        if gene_feature is None:
            group_key = ("feature", block.feature_index)
            group_feature = block
        else:
            group_key = ("gene", gene_feature.feature_index)
            group_feature = gene_feature
        groups.setdefault(group_key, []).append(block)
        group_features.setdefault(group_key, group_feature)

    sorted_group_keys = sorted(
        group_features,
        key=lambda group_key: genbank_feature_sort_key(group_features[group_key]),
    )
    sorted_blocks = [*source_blocks]
    for group_key in sorted_group_keys:
        sorted_blocks.extend(
            sorted(
                groups[group_key],
                key=lambda block: (
                    FEATURE_SORT_ORDER.get(block.key, 3),
                    genbank_feature_sort_key(block),
                ),
            )
        )

    original_feature_lines = record[features_index + 1 : origin_index]
    sorted_feature_lines = [
        line for block in sorted_blocks for line in block.lines
    ]
    changed = original_feature_lines != sorted_feature_lines
    if not changed:
        return record, False, len(blocks), len(source_blocks)

    sorted_record = [
        *record[: features_index + 1],
        *sorted_feature_lines,
        *record[origin_index:],
    ]
    return sorted_record, True, len(blocks), len(source_blocks)


def sort_genbank_features_by_location(path: Path):
    lines = path.read_text().splitlines(keepends=True)
    records = split_genbank_records(lines)
    sorted_records = []
    changed_records = 0
    feature_count = 0
    source_feature_count = 0

    for record in records:
        sorted_record, changed, record_feature_count, record_source_count = (
            sort_genbank_record_features_by_location(record)
        )
        sorted_records.append(sorted_record)
        changed_records += int(changed)
        feature_count += record_feature_count
        source_feature_count += record_source_count

    sorted_lines = [line for record in sorted_records for line in record]
    changed = lines != sorted_lines
    if changed:
        path.write_text("".join(sorted_lines))

    return {
        "feature_sort_record_count": len(records),
        "feature_sort_changed_record_count": changed_records,
        "feature_sort_feature_count": feature_count,
        "feature_sort_source_feature_count": source_feature_count,
        "feature_sort_changed": changed,
    }


def topology_from_fasta_header(header: str):
    match = FASTA_CIRCULAR_PATTERN.search(header)
    if match is None:
        return None
    return "circular" if match.group(1).lower() == "true" else "linear"


def normalize_genbank_topology(topology):
    if topology is None:
        return None
    normalized = str(topology).strip().lower()
    if normalized not in GENBANK_TOPOLOGIES:
        raise ValueError(
            f"Unsupported GenBank topology '{topology}'. Expected 'circular' or 'linear'."
        )
    return normalized


def extract_locus_metadata(line: str, sequence_length=None):
    if not line.startswith("LOCUS"):
        raise ValueError(f"Expected a GenBank LOCUS line, got: {line!r}")

    length_match = re.search(r"(\d+)\s+bp\b", line)
    if length_match:
        parsed_length = int(length_match.group(1))
        name = line[len("LOCUS") : length_match.start(1)].strip()
        suffix = line[length_match.end() :]
    else:
        parsed_length = None
        name = line[len("LOCUS") :].strip()
        suffix = ""

    length = sequence_length if sequence_length is not None else parsed_length
    if sequence_length is not None:
        length_text = str(sequence_length)
        if name.endswith(length_text):
            name = name[: -len(length_text)].rstrip()

    molecule_match = GENBANK_MOLECULE_PATTERN.search(suffix)
    molecule = molecule_match.group(1) if molecule_match else "DNA"

    topology_match = re.search(r"circular|linear", suffix, flags=re.I)
    topology = topology_match.group(0).lower() if topology_match else None

    date_match = GENBANK_DATE_PATTERN.search(suffix)
    date = date_match.group(0) if date_match else None

    division = None
    if date_match:
        before_date = suffix[: date_match.start()]
        division_match = re.search(
            r"(?:circular|linear)?\s*([A-Z]{3})\s*$",
            before_date,
            flags=re.I,
        )
        if division_match:
            division = division_match.group(1).upper()
    if division is None:
        division = "PLN"

    return {
        "name": name,
        "length": length,
        "molecule": molecule,
        "topology": topology,
        "division": division,
        "date": date,
    }


def format_locus_line(
    original_line: str,
    *,
    topology=None,
    locus_name=None,
    sequence_length=None,
):
    metadata = extract_locus_metadata(original_line, sequence_length=sequence_length)
    topology = normalize_genbank_topology(topology) or metadata["topology"]
    name = locus_name or metadata["name"]
    length = sequence_length if sequence_length is not None else metadata["length"]
    if length is None:
        raise ValueError(f"Could not determine LOCUS sequence length from: {original_line!r}")

    if topology is None:
        topology_text = ""
    else:
        topology_text = f"{topology:<8} "
    trailing = " ".join(item for item in (metadata["division"], metadata["date"]) if item)
    return (
        f"LOCUS       {name:<24} {length:>11} bp    "
        f"{metadata['molecule']:<6} {topology_text}{trailing}"
    ).rstrip()


def replace_locus_topology(line: str, topology: str):
    return format_locus_line(line, topology=topology)


def curate_genbank_locus(
    path: Path,
    *,
    topology=None,
    locus_name=None,
    sequence_length=None,
):
    topology = normalize_genbank_topology(topology)
    if topology is None and locus_name is None and sequence_length is None:
        return {
            "locus_name_before": None,
            "locus_name_after": None,
            "locus_name_changed": False,
            "locus_topology_before": None,
            "locus_topology_after": None,
            "locus_topology_changed": False,
            "locus_topology_skipped_reason": "no input FASTA circular flag provided",
            "locus_line_before": None,
            "locus_line_after": None,
            "locus_line_changed": False,
        }

    lines = path.read_text().splitlines(keepends=True)
    locus_index = None
    for index, line in enumerate(lines):
        if line.startswith("LOCUS"):
            locus_index = index
            break
    if locus_index is None:
        return {
            "locus_name_before": None,
            "locus_name_after": locus_name,
            "locus_name_changed": False,
            "locus_topology_before": None,
            "locus_topology_after": topology,
            "locus_topology_changed": False,
            "locus_topology_skipped_reason": "no LOCUS line found",
            "locus_line_before": None,
            "locus_line_after": None,
            "locus_line_changed": False,
        }

    line = lines[locus_index]
    newline = "\n" if line.endswith("\n") else ""
    line_without_newline = line.rstrip("\n")
    before_metadata = extract_locus_metadata(
        line_without_newline,
        sequence_length=sequence_length,
    )
    before_name = before_metadata["name"]
    before = before_metadata["topology"]
    after_name = locus_name or before_name
    updated_line_without_newline = format_locus_line(
        line_without_newline,
        topology=topology,
        locus_name=locus_name,
        sequence_length=sequence_length,
    )
    updated_line = updated_line_without_newline + newline
    changed = updated_line != line
    if changed:
        lines[locus_index] = updated_line
        path.write_text("".join(lines))

    return {
        "locus_name_before": before_name,
        "locus_name_after": after_name,
        "locus_name_changed": before_name != after_name,
        "locus_topology_before": before,
        "locus_topology_after": topology or before,
        "locus_topology_changed": before != (topology or before),
        "locus_topology_skipped_reason": None,
        "locus_line_before": line_without_newline,
        "locus_line_after": updated_line_without_newline,
        "locus_line_changed": changed,
    }


def curate_genbank_locus_topology(path: Path, topology):
    return curate_genbank_locus(path, topology=topology)


def curate_source_organism(path: Path, assembly_name: str):
    organism_name = organism_name_from_assembly_name(assembly_name)
    organism_line = f'{QUALIFIER_INDENT}/organism="{organism_name}"\n'
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    previous_organism = None
    organism_index = None
    for index in range(source_index + 1, next_feature_index):
        stripped = lines[index].strip()
        if stripped.startswith("/organism="):
            organism_index = index
            previous_organism = stripped.removeprefix("/organism=").strip('"')
            break

    if organism_index is None:
        lines.insert(source_index + 1, organism_line)
    elif lines[organism_index] != organism_line:
        lines[organism_index] = organism_line

    changed = previous_organism != organism_name
    if changed:
        path.write_text("".join(lines))

    return {
        "source_organism_before": previous_organism,
        "source_organism_after": organism_name,
        "source_organism_changed": changed,
    }


def curate_source_organelle(path: Path, organelle: str | None = None):
    if organelle is None:
        return {
            "source_organelle_before": None,
            "source_organelle_after": None,
            "source_organelle_changed": False,
            "source_organelle_skipped_reason": "no organelle qualifier requested",
        }
    if any(character in organelle for character in '"\r\n'):
        raise ValueError(f"Invalid source organelle qualifier: {organelle!r}")

    organelle_line = f'{QUALIFIER_INDENT}/organelle="{organelle}"\n'
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    previous_organelle = None
    organelle_index = None
    organism_index = None
    for index in range(source_index + 1, next_feature_index):
        stripped = lines[index].strip()
        if stripped.startswith("/organism="):
            organism_index = index
        elif stripped.startswith("/organelle="):
            organelle_index = index
            previous_organelle = stripped.removeprefix("/organelle=").strip('"')
            break

    if organelle_index is None:
        insert_index = organism_index + 1 if organism_index is not None else source_index + 1
        lines.insert(insert_index, organelle_line)
    elif lines[organelle_index] != organelle_line:
        lines[organelle_index] = organelle_line

    changed = previous_organelle != organelle
    if changed:
        path.write_text("".join(lines))

    return {
        "source_organelle_before": previous_organelle,
        "source_organelle_after": organelle,
        "source_organelle_changed": changed,
        "source_organelle_skipped_reason": None,
    }


def curate_source_annotation_method(path: Path, annotation_method: str | None = None):
    if annotation_method is None:
        return {
            "source_annotation_method_before": None,
            "source_annotation_method_after": None,
            "source_annotation_method_changed": False,
            "source_annotation_method_skipped_reason": (
                "no annotation method note requested"
            ),
        }
    if any(character in annotation_method for character in '"\r\n'):
        raise ValueError(f"Invalid annotation method note: {annotation_method!r}")

    note_value = f"Annotation Method :: {annotation_method}"
    note_line = f'{QUALIFIER_INDENT}/note="{note_value}"\n'
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    previous_note = None
    annotation_note_index = None
    last_source_qualifier_index = source_index
    for index in range(source_index + 1, next_feature_index):
        stripped = lines[index].strip()
        if stripped.startswith("/"):
            last_source_qualifier_index = index
        if not stripped.startswith("/note="):
            continue
        note = stripped.removeprefix("/note=").strip('"')
        if note.startswith("Annotation Method ::"):
            annotation_note_index = index
            previous_note = note
            break

    if annotation_note_index is None:
        lines.insert(last_source_qualifier_index + 1, note_line)
    elif lines[annotation_note_index] != note_line:
        lines[annotation_note_index] = note_line

    changed = previous_note != note_value
    if changed:
        path.write_text("".join(lines))

    return {
        "source_annotation_method_before": previous_note,
        "source_annotation_method_after": note_value,
        "source_annotation_method_changed": changed,
        "source_annotation_method_skipped_reason": None,
    }


def normalize_taxid(taxid):
    if taxid is None:
        return None
    if isinstance(taxid, bool):
        raise ValueError("taxid must not be a boolean.")
    value = str(taxid).strip()
    if not value:
        return None
    if not value.isdigit() or int(value) < 1:
        raise ValueError("taxid must be a positive integer NCBI Taxonomy ID.")
    return value


def curate_source_taxon_db_xref(path: Path, taxid=None, remove_without_taxid=False):
    taxid = normalize_taxid(taxid)
    desired_xref = f"taxon:{taxid}" if taxid else None
    previous_xrefs = []
    lines = path.read_text().splitlines(keepends=True)

    source_index = None
    for index, line in enumerate(lines):
        if line.startswith(f"{FEATURE_INDENT}source"):
            source_index = index
            break
    if source_index is None:
        raise ValueError(f"No source feature was found in {path}.")

    next_feature_index = len(lines)
    for index in range(source_index + 1, len(lines)):
        line = lines[index]
        if line.startswith(FEATURE_INDENT) and not line.startswith(QUALIFIER_INDENT):
            next_feature_index = index
            break

    curated_lines = []
    organism_curated_index = None
    first_taxon_curated_index = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        in_source = source_index < index < next_feature_index
        if in_source and stripped.startswith('/db_xref="taxon:'):
            if first_taxon_curated_index is None:
                first_taxon_curated_index = len(curated_lines)
            previous_xrefs.append(stripped.removeprefix("/db_xref=").strip('"'))
            continue
        if in_source and stripped.startswith("/organism="):
            organism_curated_index = len(curated_lines)
        curated_lines.append(line)

    if desired_xref is not None:
        xref_line = f'{QUALIFIER_INDENT}/db_xref="{desired_xref}"\n'
        if first_taxon_curated_index is not None:
            insert_index = first_taxon_curated_index
        elif organism_curated_index is not None:
            insert_index = organism_curated_index + 1
        else:
            insert_index = source_index + 1
        curated_lines.insert(insert_index, xref_line)
        current_xrefs = [desired_xref]
    elif previous_xrefs and not remove_without_taxid:
        insert_index = first_taxon_curated_index
        for xref in previous_xrefs:
            curated_lines.insert(
                insert_index,
                f'{QUALIFIER_INDENT}/db_xref="{xref}"\n',
            )
            insert_index += 1
        current_xrefs = previous_xrefs
    else:
        current_xrefs = []

    changed = previous_xrefs != current_xrefs
    if changed:
        path.write_text("".join(curated_lines))

    return {
        "source_taxon_db_xrefs_before": previous_xrefs,
        "source_taxon_db_xrefs_after": current_xrefs,
        "source_taxon_db_xref_changed": changed,
    }


def curate_species_placeholders(path: Path, assembly_name: str):
    organism_name = organism_name_from_assembly_name(assembly_name)
    lines = path.read_text().splitlines(keepends=True)
    replacements = 0
    curated_lines = []
    for line in lines:
        replacements += line.count(SPECIES_PLACEHOLDER)
        curated_lines.append(line.replace(SPECIES_PLACEHOLDER, organism_name))

    changed = replacements > 0
    if changed:
        path.write_text("".join(curated_lines))

    return {
        "species_placeholder": SPECIES_PLACEHOLDER,
        "species_placeholder_replacement": organism_name,
        "species_placeholder_replacements": replacements,
        "species_placeholder_changed": changed,
    }


def curate_genbank_source_metadata(
    path: Path,
    assembly_name: str,
    taxid=None,
    organelle=None,
    annotation_method=None,
):
    normalized_taxid = normalize_taxid(taxid)
    source_organism = curate_source_organism(path, assembly_name)
    source_organelle = curate_source_organelle(path, organelle)
    source_annotation_method = curate_source_annotation_method(
        path,
        annotation_method,
    )
    source_taxon_db_xref = curate_source_taxon_db_xref(
        path,
        taxid=normalized_taxid,
        remove_without_taxid=source_organism["source_organism_changed"],
    )
    species_placeholders = curate_species_placeholders(path, assembly_name)
    return {
        "source_organism": source_organism,
        "source_organelle": source_organelle,
        "source_annotation_method": source_annotation_method,
        "source_taxon_db_xref": source_taxon_db_xref,
        "species_placeholders": species_placeholders,
    }


def append_post_curation_summary(lines, post_curation):
    core_sections = post_curation.get("core_sections")
    if core_sections:
        sections = ", ".join(core_sections["core_sections"])
        record_count = core_sections["core_sections_record_count"]
        if core_sections["core_sections_changed"]:
            lines.append(
                f"- GenBank file: trimmed {record_count} record(s) to {sections}"
            )
        else:
            lines.append(
                f"- GenBank file: already contained only {sections} section(s)"
            )

    circular_origin_rotation = post_curation.get("circular_origin_rotation")
    if circular_origin_rotation:
        skipped_reason = circular_origin_rotation.get("skipped_reason")
        if circular_origin_rotation.get("rotation_applied"):
            lines.append(
                "- circular origin rotation: moved FASTA origin to position "
                f"{circular_origin_rotation['new_origin_position']} "
                f"(0-based offset {circular_origin_rotation['rotation_offset']}) "
                "inside an internal feature-free gap "
                f"{circular_origin_rotation['selected_gap_start']}.."
                f"{circular_origin_rotation['selected_gap_end']} "
                f"({circular_origin_rotation['selected_gap_length']} bp)"
            )
        elif skipped_reason:
            lines.append(f"- circular origin rotation: not applied ({skipped_reason})")
        else:
            lines.append("- circular origin rotation: not applied")

    locus_topology = post_curation.get("locus_topology")
    if locus_topology:
        before_name = locus_topology.get("locus_name_before")
        after_name = locus_topology.get("locus_name_after")
        if before_name is not None and after_name is not None:
            if locus_topology.get("locus_name_changed"):
                lines.append(
                    f'- LOCUS name: changed from "{before_name}" to "{after_name}"'
                )
            else:
                lines.append(f'- LOCUS name: already "{after_name}"')

        before = locus_topology["locus_topology_before"]
        after = locus_topology["locus_topology_after"]
        skipped_reason = locus_topology["locus_topology_skipped_reason"]
        if skipped_reason:
            lines.append(f"- LOCUS topology: skipped ({skipped_reason})")
        elif before is None:
            lines.append(f'- LOCUS topology: added "{after}"')
        elif locus_topology["locus_topology_changed"]:
            lines.append(f'- LOCUS topology: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- LOCUS topology: already "{after}"')

    origin_wrap = post_curation.get("origin_wrapping_locations")
    if origin_wrap:
        record_count = origin_wrap["origin_wrap_record_count"]
        location_count = origin_wrap["origin_wrap_location_count"]
        if origin_wrap["origin_wrap_changed"]:
            changed_records = origin_wrap["origin_wrap_changed_record_count"]
            lines.append(
                "- feature locations: normalized "
                f"{location_count} circular origin-spanning location(s) "
                f"across {changed_records}/{record_count} record(s)"
            )
        else:
            lines.append(
                "- feature locations: no circular origin-spanning location(s) "
                f"needed normalization across {record_count} record(s)"
            )

    location_wrap = post_curation.get("feature_location_wrapping")
    if location_wrap:
        record_count = location_wrap["feature_location_wrap_record_count"]
        feature_count = location_wrap["feature_location_wrap_changed_feature_count"]
        if location_wrap["feature_location_wrap_changed"]:
            changed_records = location_wrap[
                "feature_location_wrap_changed_record_count"
            ]
            lines.append(
                "- feature location wrapping: normalized "
                f"{feature_count} feature location(s) across "
                f"{changed_records}/{record_count} record(s)"
            )
        else:
            lines.append(
                "- feature location wrapping: already standard across "
                f"{record_count} record(s)"
            )

    feature_sort = post_curation.get("feature_sort")
    if feature_sort:
        record_count = feature_sort["feature_sort_record_count"]
        feature_count = feature_sort["feature_sort_feature_count"]
        source_count = feature_sort["feature_sort_source_feature_count"]
        if feature_sort["feature_sort_changed"]:
            changed_records = feature_sort["feature_sort_changed_record_count"]
            lines.append(
                "- feature table order: sorted "
                f"{feature_count} feature(s) across {changed_records}/{record_count} "
                "record(s) by genomic coordinate "
                f"(kept {source_count} source feature(s) first)"
            )
        else:
            lines.append(
                "- feature table order: already sorted "
                f"by genomic coordinate across {record_count} record(s)"
            )

    pga_cds_qc = post_curation.get("pga_v2_cds_qc")
    if pga_cds_qc:
        row_count = pga_cds_qc["pga_v2_cds_qc_row_count"]
        alignment_count = pga_cds_qc[
            "pga_v2_cds_qc_alignment_candidate_count"
        ]
        fix_count = pga_cds_qc["pga_v2_cds_qc_candidate_fix_count"]
        lines.append(
            "- PGA v2 CDS length QC: checked "
            f"{row_count} CDS feature(s); {alignment_count} short/warning "
            f"candidate(s) were aligned to reference CDS; found {fix_count} "
            "read-support-testable frameshift fix candidate(s)"
        )

    pga_sequence_fix = post_curation.get("pga_v2_sequence_frameshift_fix")
    if pga_sequence_fix:
        if pga_sequence_fix["enabled"]:
            applied_count = pga_sequence_fix["applied_fix_count"]
            if pga_sequence_fix["rerun_applied"]:
                lines.append(
                    "- PGA v2 chloroplast sequence frameshift fixes: applied "
                    f"{applied_count} high-confidence read-supported FASTA "
                    "sequence fix(es) and reran PGA v2"
                )
            else:
                lines.append(
                    "- PGA v2 chloroplast sequence frameshift fixes: enabled, "
                    "but no high-confidence read-supported FASTA sequence fixes "
                    "were applied"
                )
        else:
            lines.append(
                "- PGA v2 chloroplast sequence frameshift fixes: disabled "
                "(reported candidates only)"
            )

    translation_validation = post_curation.get("cds_auto_translation_validation")
    if translation_validation:
        checked_count = translation_validation[
            "cds_auto_translation_checked_cds_count"
        ]
        strict_table = translation_validation["cds_auto_translation_strict_table"]
        comparison_table = translation_validation[
            "cds_auto_translation_comparison_table"
        ]
        strict_pass_count = translation_validation[
            "cds_auto_translation_strict_pass_count"
        ]
        match_count = translation_validation[
            "cds_auto_translation_table_compare_match_count"
        ]
        skipped_exception_count = translation_validation[
            "cds_auto_translation_skipped_exception_count"
        ]
        skipped_partial_count = translation_validation[
            "cds_auto_translation_skipped_partial_count"
        ]
        skipped_pseudo_count = translation_validation[
            "cds_auto_translation_skipped_pseudo_count"
        ]
        blocking_count = translation_validation[
            "cds_auto_translation_blocking_issue_count"
        ]
        warning_count = translation_validation["cds_auto_translation_warning_count"]
        if translation_validation["cds_auto_translation_passed"]:
            lines.append(
                "- CDS auto-translation QC: "
                f"{strict_pass_count}/{checked_count} non-exception complete CDS "
                f"passed table {strict_table} validation; "
                f"{match_count}/{checked_count} matched between table "
                f"{strict_table} and table {comparison_table}"
            )
        else:
            lines.append(
                "- CDS auto-translation QC: found "
                f"{blocking_count} blocking issue(s) and {warning_count} "
                f"warning(s) among {checked_count} checked CDS"
            )
        if skipped_exception_count or skipped_partial_count or skipped_pseudo_count:
            lines.append(
                "- CDS auto-translation QC skipped: "
                f"{skipped_exception_count} CDS with /exception, "
                f"{skipped_partial_count} partial CDS, "
                f"{skipped_pseudo_count} pseudo/pseudogene CDS"
            )

    trans_splicing = post_curation.get("pmga_trans_splicing")
    if trans_splicing:
        record_count = trans_splicing["pmga_trans_splicing_record_count"]
        feature_keys = trans_splicing.get("pmga_trans_splicing_feature_keys", ["CDS"])
        feature_text = "/".join(feature_keys)
        changed_features = trans_splicing["pmga_trans_splicing_changed_feature_count"]
        removed_count = trans_splicing[
            "pmga_trans_splicing_removed_exception_count"
        ]
        added_count = trans_splicing["pmga_trans_splicing_added_qualifier_count"]
        if trans_splicing["pmga_trans_splicing_changed"]:
            changed_records = trans_splicing[
                "pmga_trans_splicing_changed_record_count"
            ]
            lines.append(
                "- PMGA trans-splicing qualifiers: converted "
                f'{removed_count} /exception="trans-splicing" qualifier(s) '
                f"to {added_count} /trans_splicing qualifier(s) on "
                f"{changed_features} {feature_text} feature(s) across "
                f"{changed_records}/{record_count} record(s)"
            )
        else:
            lines.append(
                "- PMGA trans-splicing qualifiers: no "
                f'/exception="trans-splicing" qualifier(s) needed normalization '
                f"across {record_count} record(s)"
            )

    qualifier_removal = post_curation.get("feature_qualifier_removal")
    if qualifier_removal:
        feature_keys = qualifier_removal["qualifier_removal_feature_keys"]
        qualifiers = qualifier_removal["qualifier_removal_qualifiers"]
        qualifier_counts = qualifier_removal["qualifier_removal_counts"]
        record_count = qualifier_removal["qualifier_removal_record_count"]
        target_feature_count = qualifier_removal[
            "qualifier_removal_target_feature_count"
        ]
        feature_text = "/".join(feature_keys)
        qualifier_text = ", ".join(f"/{qualifier}" for qualifier in qualifiers)
        counted_qualifier_text = ", ".join(
            f"/{qualifier} ({qualifier_counts.get(qualifier, 0)})"
            for qualifier in qualifiers
        )
        if qualifier_removal["qualifier_removal_changed"]:
            changed_records = qualifier_removal[
                "qualifier_removal_changed_record_count"
            ]
            changed_features = qualifier_removal[
                "qualifier_removal_changed_feature_count"
            ]
            lines.append(
                f"- {feature_text} qualifiers: removed {counted_qualifier_text} "
                f"from {changed_features}/{target_feature_count} "
                f"{feature_text} feature(s) across {changed_records}/{record_count} "
                "record(s)"
            )
            if "transl_table" in qualifiers:
                lines.append(
                    "- CDS /transl_table: omitted to use the INSDC default "
                    "standard genetic code (table 1), equivalent to explicitly "
                    "setting /transl_table=1"
                )
        else:
            lines.append(
                f"- {feature_text} qualifiers: no {qualifier_text} qualifier(s) "
                f"needed removal across {target_feature_count} feature(s)"
            )
            if "transl_table" in qualifiers:
                lines.append(
                    "- CDS /transl_table: omitted qualifiers remain interpreted "
                    "as the INSDC default standard genetic code (table 1)"
                )

    source_organism = post_curation.get("source_organism")
    if source_organism:
        before = source_organism["source_organism_before"]
        after = source_organism["source_organism_after"]
        if before is None:
            lines.append(f'- source /organism: added "{after}"')
        elif source_organism["source_organism_changed"]:
            lines.append(f'- source /organism: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- source /organism: already "{after}"')

    source_organelle = post_curation.get("source_organelle")
    if source_organelle and not source_organelle["source_organelle_skipped_reason"]:
        before = source_organelle["source_organelle_before"]
        after = source_organelle["source_organelle_after"]
        if before is None:
            lines.append(f'- source /organelle: added "{after}"')
        elif source_organelle["source_organelle_changed"]:
            lines.append(f'- source /organelle: changed from "{before}" to "{after}"')
        else:
            lines.append(f'- source /organelle: already "{after}"')

    source_annotation_method = post_curation.get("source_annotation_method")
    if (
        source_annotation_method
        and not source_annotation_method["source_annotation_method_skipped_reason"]
    ):
        before = source_annotation_method["source_annotation_method_before"]
        after = source_annotation_method["source_annotation_method_after"]
        if before is None:
            lines.append(f'- source annotation method /note: added "{after}"')
        elif source_annotation_method["source_annotation_method_changed"]:
            lines.append(
                f'- source annotation method /note: changed from "{before}" '
                f'to "{after}"'
            )
        else:
            lines.append(f'- source annotation method /note: already "{after}"')

    source_taxon_db_xref = post_curation.get("source_taxon_db_xref")
    if source_taxon_db_xref and source_taxon_db_xref["source_taxon_db_xref_changed"]:
        before = source_taxon_db_xref["source_taxon_db_xrefs_before"]
        after = source_taxon_db_xref["source_taxon_db_xrefs_after"]
        before_text = ", ".join(before) if before else "none"
        after_text = ", ".join(after) if after else "none"
        lines.append(
            f"- source taxon db_xref: changed from {before_text} to {after_text}"
        )

    species_placeholders = post_curation.get("species_placeholders")
    if species_placeholders and species_placeholders["species_placeholder_changed"]:
        placeholder = species_placeholders["species_placeholder"]
        replacement = species_placeholders["species_placeholder_replacement"]
        count = species_placeholders["species_placeholder_replacements"]
        lines.append(
            f'- species placeholders: replaced {count} "{placeholder}" '
            f'occurrence(s) with "{replacement}"'
        )


def format_post_curation_record(post_curation, annotation_path: Path | None = None):
    if annotation_path is None:
        manual_note_line = (
            "This file is regenerated by the workflow. If you edit the GenBank "
            "annotation after the organelle annotation step finishes, record those "
            "edits below."
        )
    else:
        manual_note_line = (
            "This file is regenerated by the workflow. If you edit "
            f"`{annotation_path.name}` after the organelle annotation step "
            "finishes, record those edits below."
        )

    lines = [
        manual_note_line,
        "",
        "# Automatic post-curation",
        "",
    ]

    record_curations = post_curation.get("records")
    if record_curations:
        lines.append(
            "- PMGA was run separately for each input FASTA record; the curated "
            "GenBank records were then concatenated into the final multi-record file."
        )
        for record_curation in record_curations:
            record_id = record_curation.get("input_record_id")
            record_annotation = record_curation.get("annotation")
            label = f"`{record_id}`" if record_id else "record"
            if record_annotation:
                label = f"{label} -> `{Path(record_annotation).name}`"
            lines.append(f"- {label}")
        lines.append("")

        for record_curation in record_curations:
            record_id = record_curation.get("input_record_id", "record")
            lines.extend([f"## {record_id}", ""])
            append_post_curation_summary(
                lines,
                record_curation.get("post_curation", {}),
            )
            lines.append("")
    else:
        append_post_curation_summary(lines, post_curation)

    lines.extend(
        [
            "",
            "# Manual post-curation",
            "",
        ]
    )
    return "\n".join(lines)


def write_post_curation_record(
    path: Path,
    post_curation,
    annotation_path: Path | None = None,
):
    updated = format_post_curation_record(post_curation, annotation_path)

    path.parent.mkdir(parents=True, exist_ok=True)
    changed = (not path.exists()) or path.read_text() != updated
    if changed:
        path.write_text(updated)

    return {
        "post_curation_record": str(path),
        "post_curation_record_changed": changed,
    }


def write_run_manifest(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
