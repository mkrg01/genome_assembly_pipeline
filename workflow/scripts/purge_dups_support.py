"""Validation and conservative extraction policy for the optional purge_dups stage."""

import argparse
import json
from pathlib import Path
import subprocess
import sys


def fasta_lengths(path):
    """Read lengths without retaining chromosome sequences in memory."""
    lengths = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or ":" in fields[0] or fields[0] in lengths:
                    raise ValueError(f"{path}: empty, duplicate, or colon-containing FASTA ID")
                name = fields[0]
                lengths[name] = 0
            elif name is None:
                raise ValueError(f"{path}: sequence before FASTA header")
            else:
                lengths[name] += len(line)
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError(f"{path}: empty assembly or empty FASTA record")
    return lengths


def read_cutoffs(text):
    try:
        values = [int(value) for value in text.split()]
    except ValueError as error:
        raise ValueError("Invalid calcuts output: expected six integer cutoffs") from error
    if (len(values) != 6 or values != sorted(values) or values[0] < 0
            or not values[0] < values[3] < values[5]):
        raise ValueError(f"Invalid calcuts output: {values}; inspect coverage and set purge_dups_cutoffs manually")
    return values


def resolve_cutoffs(stat, output, metadata, manual=None):
    command = ["calcuts"]
    if manual is not None:
        low, mid, high = manual
        if not 0 <= low < mid < high:
            raise ValueError("Manual cutoffs must satisfy 0 <= low < mid < high")
        command.extend(["-l", str(low), "-m", str(mid), "-u", str(high)])
    command.append(str(stat))
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    print(result.stderr, file=sys.stderr, end="")
    values = read_cutoffs(result.stdout)
    warnings = [line for line in result.stderr.splitlines() if "[W::" in line]
    Path(output).write_text("\t".join(map(str, values)) + "\n")
    Path(metadata).write_text(json.dumps({
        "mode": "manual" if manual is not None else "auto",
        "cutoffs": values, "warnings": warnings, "command": command,
    }, indent=2) + "\n")


def filter_bed(assembly, bed, output):
    """Keep haplotypic calls; depth-only JUNK/HIGHCOV/REPEAT calls stay in raw BED."""
    lengths = fasta_lengths(assembly)
    records = []
    with Path(bed).open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"{bed}:{number}: invalid BED record")
            name, start, end, kind = fields[:4]
            start, end = int(start), int(end)
            if name not in lengths or not 0 <= start < end <= lengths[name]:
                raise ValueError(f"{bed}:{number}: coordinates outside input assembly")
            if kind in {"HAPLOTIG", "OVLP"}:
                records.append((name, start, end, kind))
    # get_seqs expects records grouped by contig and sorted by position.
    with Path(output).open("w") as handle:
        for record in sorted(set(records)):
            handle.write("\t".join(map(str, record)) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    size = sub.add_parser("size")
    size.add_argument("assembly")
    cutoffs = sub.add_parser("cutoffs")
    for name in ("stat", "output", "metadata"):
        cutoffs.add_argument("--" + name, required=True)
    cutoffs.add_argument("--manual", nargs=3, type=int)
    bed = sub.add_parser("filter-bed")
    for name in ("assembly", "bed", "output"):
        bed.add_argument("--" + name, required=True)
    args = vars(parser.parse_args())
    action = args.pop("action")
    try:
        if action == "size":
            print(sum(fasta_lengths(args["assembly"]).values()))
        elif action == "cutoffs":
            resolve_cutoffs(**args)
        elif action == "filter-bed":
            filter_bed(**args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
