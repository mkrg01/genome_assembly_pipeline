import argparse
import csv
from pathlib import Path


def normalize_int(value):
    if value is None:
        raise ValueError("sum_len is empty")
    normalized = value.replace(",", "").strip()
    if not normalized.isdigit():
        raise ValueError(f"sum_len must be a positive integer, got {value!r}")
    parsed = int(normalized)
    if parsed <= 0:
        raise ValueError(f"sum_len must be greater than zero, got {value!r}")
    return parsed


def extract_sum_len(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty or missing a header")
        if "sum_len" not in reader.fieldnames:
            raise ValueError(f"{path} is missing the required sum_len column")
        try:
            row = next(reader)
        except StopIteration as error:
            raise ValueError(f"{path} has a header but no data rows") from error
    return normalize_int(row["sum_len"])


def main():
    parser = argparse.ArgumentParser(
        description="Extract the sum_len value from a seqkit stats TSV file."
    )
    parser.add_argument("seqkit_stats", type=Path)
    args = parser.parse_args()
    try:
        print(extract_sum_len(args.seqkit_stats))
    except ValueError as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
