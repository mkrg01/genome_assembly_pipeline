#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path


PRIMARY_ORGANELLE_EXCLUDE_FLAGS = 4 + 256 + 512 + 2048


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract primary alignments on one organelle's contigs from a "
            "nuclear+organelle coordinate-sorted BAM."
        )
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--organelle", required=True)
    parser.add_argument("--input-bam", type=Path, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    return parser.parse_args()


def organelle_record_ids(manifest_path: Path, organelle: str):
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest.get("organelles", []):
        if entry.get("name") == organelle:
            record_ids = entry.get("record_ids", [])
            if not record_ids:
                raise RuntimeError(
                    f"Manifest has no record IDs for organelle {organelle!r}."
                )
            return record_ids
    names = ", ".join(
        sorted(entry.get("name", "?") for entry in manifest.get("organelles", []))
    )
    raise RuntimeError(
        f"Organelle {organelle!r} was not found in {manifest_path}. "
        f"Available organelles: {names or 'none'}."
    )


def main():
    args = parse_args()
    record_ids = organelle_record_ids(args.manifest, args.organelle)
    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "samtools",
            "view",
            "-b",
            "-F",
            str(PRIMARY_ORGANELLE_EXCLUDE_FLAGS),
            "-o",
            str(args.output_bam),
            str(args.input_bam),
            *record_ids,
        ],
        check=True,
    )
    subprocess.run(["samtools", "index", "-f", str(args.output_bam)], check=True)


if __name__ == "__main__":
    main()
