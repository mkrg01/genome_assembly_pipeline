#!/usr/bin/env python3
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


COPY_BUFFER_SIZE = 1024 * 1024


@dataclass(frozen=True)
class FastaRecordIndex:
    original_id: str
    original_header: str
    sequence_start: int
    record_end: int
    length: int


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Sort assembly FASTA records by decreasing length and rename them "
            "using a sequential prefix."
        )
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--source-stage", required=True)
    return parser.parse_args()


def _decode_header(raw_header: bytes, path: Path, line_number: int):
    try:
        header = raw_header.decode("utf-8").strip()
    except UnicodeDecodeError as error:
        raise ValueError(
            f"FASTA header in {path} at line {line_number} is not valid UTF-8."
        ) from error
    if not header:
        raise ValueError(f"Empty FASTA header found in {path} at line {line_number}.")
    return header


def index_fasta(path: Path):
    records = []
    seen_ids = set()
    current = None
    offset = 0

    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            line_start = offset
            offset += len(line)

            if line.startswith(b">"):
                if current is not None:
                    records.append(
                        FastaRecordIndex(
                            original_id=current["original_id"],
                            original_header=current["original_header"],
                            sequence_start=current["sequence_start"],
                            record_end=line_start,
                            length=current["length"],
                        )
                    )

                header = _decode_header(line[1:].strip(), path, line_number)
                original_id = header.split(maxsplit=1)[0]
                if original_id in seen_ids:
                    raise ValueError(
                        f"Duplicate FASTA record ID {original_id!r} found in {path}."
                    )
                seen_ids.add(original_id)
                current = {
                    "original_id": original_id,
                    "original_header": header,
                    "sequence_start": offset,
                    "length": 0,
                }
            elif current is None:
                if line.strip():
                    raise ValueError(
                        f"Sequence data found before the first FASTA header in {path} "
                        f"at line {line_number}."
                    )
            else:
                current["length"] += len(b"".join(line.split()))

    if current is not None:
        records.append(
            FastaRecordIndex(
                original_id=current["original_id"],
                original_header=current["original_header"],
                sequence_start=current["sequence_start"],
                record_end=offset,
                length=current["length"],
            )
        )

    if not records:
        raise ValueError(f"No FASTA records found in {path}.")
    empty_ids = [record.original_id for record in records if record.length == 0]
    if empty_ids:
        raise ValueError(
            f"Empty FASTA sequence found in {path} for record {empty_ids[0]!r}."
        )
    return records


def validate_prefix(prefix: str):
    if not prefix or any(character.isspace() for character in prefix):
        raise ValueError("The FASTA record prefix must be non-empty and contain no whitespace.")


def _copy_sequence(input_handle, output_handle, record: FastaRecordIndex):
    input_handle.seek(record.sequence_start)
    remaining = record.record_end - record.sequence_start
    last_byte = b""
    while remaining:
        chunk = input_handle.read(min(COPY_BUFFER_SIZE, remaining))
        if not chunk:
            raise OSError("Unexpected end of FASTA while copying a sequence record.")
        output_handle.write(chunk)
        last_byte = chunk[-1:]
        remaining -= len(chunk)
    if last_byte not in {b"\n", b"\r"}:
        output_handle.write(b"\n")


def rename_assembly(input_path, output_path, mapping_path, prefix, source_stage):
    validate_prefix(prefix)
    records = sorted(
        index_fasta(input_path),
        key=lambda record: (-record.length, record.original_id),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    with (
        input_path.open("rb") as input_handle,
        output_path.open("wb") as output_handle,
        mapping_path.open("w", newline="") as mapping_handle,
    ):
        mapping_writer = csv.writer(mapping_handle, delimiter="\t", lineterminator="\n")
        mapping_writer.writerow(
            [
                "rank",
                "new_name",
                "original_id",
                "original_header",
                "length",
                "source_stage",
            ]
        )
        for rank, record in enumerate(records, start=1):
            new_name = f"{prefix}{rank}"
            output_handle.write(f">{new_name}\n".encode("utf-8"))
            _copy_sequence(input_handle, output_handle, record)
            mapping_writer.writerow(
                [
                    rank,
                    new_name,
                    record.original_id,
                    record.original_header,
                    record.length,
                    source_stage,
                ]
            )


def main():
    args = parse_args()
    rename_assembly(
        args.input,
        args.output,
        args.mapping,
        args.prefix,
        args.source_stage,
    )


if __name__ == "__main__":
    main()
