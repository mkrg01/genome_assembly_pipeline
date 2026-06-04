import argparse
import shutil
import subprocess
from pathlib import Path

from organelle_annotation_utils import (
    copy_first_genbank,
    curate_genbank_species_metadata,
    load_taxonomy_lineage_record,
    write_post_curation_record,
    write_run_manifest,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Run PGA v2.0 on an Oatk plastid FASTA.")
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--post-curation", type=Path, required=True)
    parser.add_argument("--assembly-name", required=True)
    parser.add_argument("--form", default="circular")
    parser.add_argument("--ir", default="1000")
    parser.add_argument("--pidentity", default="40")
    parser.add_argument("--link", default="Y")
    parser.add_argument("--redundancy", default="N")
    parser.add_argument("--qcoverage", default="0.5,2.0")
    parser.add_argument("--warning", default="warning")
    parser.add_argument("--taxid")
    parser.add_argument("--taxonomy-lineage")
    return parser.parse_args()


def recreate_dir(path):
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def main():
    args = parse_args()
    script = args.script.resolve()
    input_fasta = args.input_fasta.resolve()
    reference_dir = args.reference_dir.resolve()
    annotation = args.annotation.resolve()
    manifest = args.manifest.resolve()
    post_curation_path = args.post_curation.resolve()

    if not reference_dir.is_dir():
        raise RuntimeError(
            f"PGA v2 reference directory does not exist: {reference_dir}. "
            "Set 'pga_v2_reference_dir' in config/config.yml."
        )

    run_root = annotation.parent / "pga_v2_work"
    target_dir = run_root / "target"
    copied_reference_dir = run_root / "reference"
    raw_output_dir = run_root / "gb"
    recreate_dir(run_root)
    target_dir.mkdir(parents=True)
    shutil.copytree(reference_dir, copied_reference_dir)
    staged_fasta = target_dir / f"{input_fasta.stem}.fa"
    shutil.copy2(input_fasta, staged_fasta)

    cmd = [
        "perl",
        str(script),
        "-r",
        str(copied_reference_dir),
        "-t",
        str(target_dir),
        "-i",
        str(args.ir),
        "-p",
        str(args.pidentity),
        "-l",
        args.link,
        "-d",
        args.redundancy,
        "-q",
        args.qcoverage,
        "-o",
        str(raw_output_dir),
        "-f",
        args.form,
        "-w",
        args.warning,
    ]
    subprocess.run(cmd, check=True, cwd=run_root)
    selected = copy_first_genbank(raw_output_dir, annotation)
    taxonomy_lineage = load_taxonomy_lineage_record(args.taxonomy_lineage)
    post_curation = curate_genbank_species_metadata(
        annotation,
        args.assembly_name,
        taxid=args.taxid,
        taxonomy_lineage=taxonomy_lineage,
    )
    post_curation_record = write_post_curation_record(
        post_curation_path,
        post_curation,
        annotation,
    )
    write_run_manifest(
        manifest,
        {
            "tool": "pga_v2",
            "script": str(script),
            "input_fasta": str(input_fasta),
            "reference_dir": str(reference_dir),
            "staged_reference_dir": str(copied_reference_dir),
            "raw_output_dir": str(raw_output_dir),
            "selected_annotation": str(selected),
            "annotation": str(annotation),
            "taxonomy_lineage": args.taxonomy_lineage or None,
            "post_curation": post_curation,
            **post_curation_record,
            "command": cmd,
        },
    )


if __name__ == "__main__":
    main()
