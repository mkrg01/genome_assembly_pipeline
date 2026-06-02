import argparse
import shutil
import subprocess
from pathlib import Path

from organelle_annotation_utils import copy_first_genbank, write_run_manifest


def parse_args():
    parser = argparse.ArgumentParser(description="Run MitoZ annotation on an Oatk mitochondrial FASTA.")
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--clade", required=True)
    parser.add_argument("--genetic-code", required=True)
    parser.add_argument("--threads", type=int, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    input_fasta = args.input_fasta.resolve()
    annotation = args.annotation.resolve()
    manifest = args.manifest.resolve()
    raw_output_dir = annotation.parent / "mitoz_raw"
    if raw_output_dir.exists():
        shutil.rmtree(raw_output_dir)
    raw_output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mitoz",
        "annotate",
        "--outprefix",
        args.prefix,
        "--thread_number",
        str(args.threads),
        "--clade",
        args.clade,
        "--genetic_code",
        str(args.genetic_code),
        "--fastafile",
        str(input_fasta),
    ]
    subprocess.run(cmd, check=True, cwd=raw_output_dir)
    selected = copy_first_genbank(raw_output_dir, annotation)
    write_run_manifest(
        manifest,
        {
            "tool": "mitoz",
            "input_fasta": str(input_fasta),
            "raw_output_dir": str(raw_output_dir),
            "selected_annotation": str(selected),
            "annotation": str(annotation),
            "clade": args.clade,
            "genetic_code": str(args.genetic_code),
            "command": cmd,
        },
    )


if __name__ == "__main__":
    main()
