import argparse
import gzip
import itertools
import json
import math
from pathlib import Path

from extract_seqkit_sum_len import extract_sum_len


AUTO_MINIMUM_KMER_COVERAGE_MULTIPLIER = 7


def normalize_configured_value(value):
    if value is None:
        raise ValueError("oatk_minimum_kmer_coverage is required")
    value = str(value).strip()
    if value.lower() == "auto":
        return "auto"
    if not value.isdigit():
        raise ValueError(
            "oatk_minimum_kmer_coverage must be a positive integer or 'auto'"
        )
    parsed = int(value)
    if parsed <= 0:
        raise ValueError(
            "oatk_minimum_kmer_coverage must be greater than zero or set to 'auto'"
        )
    return str(parsed)


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def first_nonempty_line(handle):
    for line in handle:
        if line.strip():
            return line
    return None


def count_fasta_bases(first_line, handle):
    total_bases = 0
    record_count = 0
    in_record = False
    for line in itertools.chain([first_line], handle):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            record_count += 1
            in_record = True
            continue
        if not in_record:
            raise ValueError("FASTA sequence data appeared before the first header")
        total_bases += len(line)
    if record_count == 0 or total_bases <= 0:
        raise ValueError("input FASTA contains no sequence bases")
    return total_bases


def count_fastq_bases(first_line, handle):
    total_bases = 0
    record_count = 0
    header = first_line
    while header:
        if not header.startswith("@"):
            raise ValueError(f"FASTQ header must start with '@', got {header!r}")
        sequence = handle.readline()
        plus = handle.readline()
        _quality = handle.readline()
        if not sequence or not plus or not _quality:
            raise ValueError("truncated FASTQ record")
        if not plus.startswith("+"):
            raise ValueError(f"FASTQ separator must start with '+', got {plus!r}")
        total_bases += len(sequence.strip())
        record_count += 1
        header = handle.readline()
    if record_count == 0 or total_bases <= 0:
        raise ValueError("input FASTQ contains no sequence bases")
    return total_bases


def count_sequence_bases(path):
    path = Path(path)
    with open_text(path) as handle:
        first_line = first_nonempty_line(handle)
        if first_line is None:
            raise ValueError(f"{path} is empty")
        if first_line.startswith(">"):
            return count_fasta_bases(first_line, handle)
        if first_line.startswith("@"):
            return count_fastq_bases(first_line, handle)
    raise ValueError(f"{path} is not recognized as FASTA or FASTQ")


def resolve_manual(configured_value):
    resolved = int(configured_value)
    return {
        "configured_value": configured_value,
        "mode": "manual",
        "resolved_minimum_kmer_coverage": resolved,
        "auto_multiplier": None,
        "hifi_read_bases": None,
        "nuclear_assembly_bases": None,
        "nuclear_coverage": None,
        "formula": "manual",
        "source_files": {},
    }


def resolve_auto(hifi_reads, nuclear_assembly_stats):
    if hifi_reads is None or nuclear_assembly_stats is None:
        raise ValueError(
            "'auto' mode requires both --hifi-reads and --nuclear-assembly-stats"
        )
    hifi_read_bases = count_sequence_bases(hifi_reads)
    nuclear_assembly_bases = extract_sum_len(nuclear_assembly_stats)
    nuclear_coverage = hifi_read_bases / nuclear_assembly_bases
    resolved = max(
        1,
        math.ceil(nuclear_coverage * AUTO_MINIMUM_KMER_COVERAGE_MULTIPLIER),
    )
    return {
        "configured_value": "auto",
        "mode": "auto",
        "resolved_minimum_kmer_coverage": resolved,
        "auto_multiplier": AUTO_MINIMUM_KMER_COVERAGE_MULTIPLIER,
        "hifi_read_bases": hifi_read_bases,
        "nuclear_assembly_bases": nuclear_assembly_bases,
        "nuclear_coverage": nuclear_coverage,
        "formula": (
            "ceil(hifi_read_bases / nuclear_assembly_bases * "
            f"{AUTO_MINIMUM_KMER_COVERAGE_MULTIPLIER})"
        ),
        "source_files": {
            "hifi_reads": str(hifi_reads),
            "nuclear_assembly_stats": str(nuclear_assembly_stats),
        },
    }


def resolve_oatk_minimum_kmer_coverage(
    configured_value,
    hifi_reads=None,
    nuclear_assembly_stats=None,
):
    configured_value = normalize_configured_value(configured_value)
    if configured_value == "auto":
        return resolve_auto(hifi_reads, nuclear_assembly_stats)
    return resolve_manual(configured_value)


def write_resolution(resolution, output, metadata):
    output = Path(output)
    metadata = Path(metadata)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(f"{resolution['resolved_minimum_kmer_coverage']}\n")
    metadata.write_text(json.dumps(resolution, indent=2, sort_keys=True) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resolve oatk_minimum_kmer_coverage to a positive integer."
    )
    parser.add_argument("--config-value", required=True)
    parser.add_argument("--hifi-reads", type=Path)
    parser.add_argument("--nuclear-assembly-stats", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        resolution = resolve_oatk_minimum_kmer_coverage(
            args.config_value,
            hifi_reads=args.hifi_reads,
            nuclear_assembly_stats=args.nuclear_assembly_stats,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error
    write_resolution(resolution, args.output, args.metadata)
    print(
        "Resolved oatk_minimum_kmer_coverage="
        f"{resolution['resolved_minimum_kmer_coverage']}"
    )


if __name__ == "__main__":
    main()
