#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a combined nuclear + organelle reference for organelle RNA "
            "editing read mapping."
        )
    )
    parser.add_argument("--nuclear", type=Path, required=True)
    parser.add_argument(
        "--organelle",
        action="append",
        default=[],
        metavar="NAME=FASTA",
        help="Organelle FASTA to append while preserving record IDs.",
    )
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def fasta_records(path: Path):
    header = None
    sequence_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
    if header is not None:
        yield header, "".join(sequence_parts)


def fasta_record_ids(path: Path):
    ids = []
    for header, _sequence in fasta_records(path):
        ids.append(header.split()[0])
    return ids


def parse_organelle_args(values):
    organelles = []
    for value in values:
        if "=" not in value:
            raise ValueError(
                f"Invalid --organelle value {value!r}. Expected NAME=FASTA."
            )
        name, path = value.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(
                f"Invalid --organelle value {value!r}. Organelle name is empty."
            )
        organelles.append((name, Path(path)))
    return organelles


def write_record(handle, header, sequence):
    handle.write(f">{header}\n")
    normalized = sequence.upper().replace("U", "T")
    for start in range(0, len(normalized), 80):
        handle.write(normalized[start : start + 80] + "\n")


def unique_nuclear_id(record_id, seen_ids):
    if record_id not in seen_ids:
        return record_id
    base = f"nuclear__{record_id}"
    candidate = base
    suffix = 2
    while candidate in seen_ids:
        candidate = f"{base}_{suffix}"
        suffix += 1
    return candidate


def main():
    args = parse_args()
    organelles = parse_organelle_args(args.organelle)

    organelle_ids = {}
    for name, path in organelles:
        ids = fasta_record_ids(path)
        for record_id in ids:
            if record_id in organelle_ids:
                raise ValueError(
                    "Duplicate organelle FASTA record ID "
                    f"{record_id!r} in {path} and {organelle_ids[record_id]!r}."
                )
            organelle_ids[record_id] = str(path)

    args.reference.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "nuclear": {
            "path": str(args.nuclear),
            "record_count": 0,
            "renamed_record_count": 0,
            "renamed_records": [],
        },
        "organelles": [],
        "reference": str(args.reference),
    }
    seen_ids = set(organelle_ids)

    with args.reference.open("w") as output:
        nuclear_seen = set()
        for header, sequence in fasta_records(args.nuclear):
            record_id = header.split()[0]
            output_id = unique_nuclear_id(record_id, seen_ids | nuclear_seen)
            nuclear_seen.add(output_id)
            manifest["nuclear"]["record_count"] += 1
            if output_id != record_id:
                manifest["nuclear"]["renamed_record_count"] += 1
                manifest["nuclear"]["renamed_records"].append(
                    {"from": record_id, "to": output_id}
                )
            rest = header[len(record_id) :].strip()
            output_header = output_id if not rest else f"{output_id} {rest}"
            write_record(output, output_header, sequence)

        for name, path in organelles:
            organelle_manifest = {
                "name": name,
                "path": str(path),
                "record_count": 0,
                "record_ids": [],
            }
            for header, sequence in fasta_records(path):
                record_id = header.split()[0]
                if record_id in nuclear_seen:
                    raise ValueError(
                        f"Organelle FASTA record ID {record_id!r} collides with "
                        "a nuclear record after nuclear renaming."
                    )
                organelle_manifest["record_count"] += 1
                organelle_manifest["record_ids"].append(record_id)
                write_record(output, header, sequence)
            manifest["organelles"].append(organelle_manifest)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
