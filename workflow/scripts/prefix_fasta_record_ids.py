#!/usr/bin/env python3
import argparse
import json
import textwrap
from pathlib import Path


LINE_WIDTH = 80


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prefix FASTA record IDs while preserving sequence and header metadata."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--label", required=True)
    return parser.parse_args()


def fasta_records(path: Path):
    header = None
    sequence_parts = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Empty FASTA header found in {path} at line {line_number}."
                    )
                sequence_parts = []
            elif header is not None:
                sequence_parts.append("".join(line.split()))
            elif line.strip():
                raise ValueError(
                    f"Sequence data found before the first FASTA header in {path} "
                    f"at line {line_number}."
                )
    if header is not None:
        yield header, "".join(sequence_parts)


def prefixed_header(header: str, prefix: str):
    parts = header.split(maxsplit=1)
    record_id = parts[0]
    rest = parts[1] if len(parts) == 2 else ""
    prefixed_id = record_id if record_id.startswith(prefix) else f"{prefix}{record_id}"
    return (
        prefixed_id if not rest else f"{prefixed_id} {rest}",
        record_id,
        prefixed_id,
    )


def write_record(handle, header: str, sequence: str):
    handle.write(f">{header}\n")
    wrapped = "\n".join(textwrap.wrap(sequence.upper().replace("U", "T"), LINE_WIDTH))
    if wrapped:
        handle.write(wrapped + "\n")


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen_ids = set()
    records = []
    with args.output.open("w") as output:
        for header, sequence in fasta_records(args.input):
            output_header, original_id, prefixed_id = prefixed_header(
                header,
                args.prefix,
            )
            if prefixed_id in seen_ids:
                raise ValueError(
                    f"Duplicate prefixed FASTA record ID {prefixed_id!r} in {args.input}."
                )
            seen_ids.add(prefixed_id)
            write_record(output, output_header, sequence)
            records.append(
                {
                    "original_id": original_id,
                    "prefixed_id": prefixed_id,
                    "changed": original_id != prefixed_id,
                    "length": len(sequence),
                }
            )

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(
        json.dumps(
            {
                "label": args.label,
                "prefix": args.prefix,
                "input": str(args.input),
                "output": str(args.output),
                "record_count": len(records),
                "records": records,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
