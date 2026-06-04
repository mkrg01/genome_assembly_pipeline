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


IMAGE_SUFFIXES = (".sif", ".simg")


def parse_args():
    parser = argparse.ArgumentParser(description="Run PMGA on an Oatk mitochondrial FASTA.")
    parser.add_argument("--pmga-bundle", type=Path, required=True)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--post-curation", type=Path, required=True)
    parser.add_argument("--db", default="1")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--assembly-name")
    parser.add_argument("--taxid")
    parser.add_argument("--taxonomy-lineage")
    return parser.parse_args()


def find_container_runtime():
    for name in ("apptainer", "singularity"):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError("Could not find 'apptainer' or 'singularity' on PATH.")


def resolve_pmga_container(bundle):
    bundle = bundle.resolve()
    if bundle.is_file() and bundle.suffix.lower() in IMAGE_SUFFIXES:
        return bundle
    if not bundle.is_dir():
        raise RuntimeError(f"PMGA bundle does not exist: {bundle}")

    images = sorted(
        path
        for suffix in IMAGE_SUFFIXES
        for path in bundle.rglob(f"*{suffix}")
        if path.is_file()
    )
    if images:
        return images[0]

    child_dirs = [child for child in bundle.iterdir() if child.is_dir()]
    for candidate in [bundle, *child_dirs]:
        if (candidate / "bin").exists() or (candidate / "usr").exists():
            return candidate
    return bundle


def main():
    args = parse_args()
    runtime = find_container_runtime()
    container = resolve_pmga_container(args.pmga_bundle)
    input_fasta = args.input_fasta.resolve()
    annotation = args.annotation.resolve()
    manifest = args.manifest.resolve()
    post_curation_path = args.post_curation.resolve()
    raw_output_dir = annotation.parent / "pmga_raw"
    if raw_output_dir.exists():
        shutil.rmtree(raw_output_dir)
    raw_output_dir.mkdir(parents=True, exist_ok=True)

    cwd = Path.cwd().resolve()
    cmd = [
        runtime,
        "exec",
        "--bind",
        f"{cwd}:{cwd}",
        str(container),
        "mgavas_m",
        "-pid",
        args.prefix,
        "--in",
        str(input_fasta),
        "--db",
        str(args.db),
        "--outdir",
        str(raw_output_dir),
    ]
    subprocess.run(cmd, check=True)
    selected = copy_first_genbank(raw_output_dir, annotation)
    assembly_name = args.assembly_name or args.prefix
    taxonomy_lineage = load_taxonomy_lineage_record(args.taxonomy_lineage)
    post_curation = curate_genbank_species_metadata(
        annotation,
        assembly_name,
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
            "tool": "pmga",
            "container_runtime": runtime,
            "pmga_bundle": str(args.pmga_bundle.resolve()),
            "pmga_container": str(container),
            "input_fasta": str(input_fasta),
            "raw_output_dir": str(raw_output_dir),
            "selected_annotation": str(selected),
            "annotation": str(annotation),
            "taxonomy_lineage": args.taxonomy_lineage or None,
            "post_curation": post_curation,
            **post_curation_record,
            "db": str(args.db),
            "command": cmd,
        },
    )


if __name__ == "__main__":
    main()
