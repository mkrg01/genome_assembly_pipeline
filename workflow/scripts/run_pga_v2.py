import argparse
import shutil
import subprocess
from pathlib import Path

from organelle_annotation_utils import (
    copy_first_genbank,
    curate_genbank_locus,
    curate_genbank_source_metadata,
    normalize_genbank_origin_wrapping_locations,
    parse_fasta_records,
    rotate_fasta_record,
    select_circular_origin_rotation_from_genbank,
    sort_genbank_features_by_location,
    topology_from_fasta_header,
    trim_genbank_to_core_sections,
    write_fasta_records,
    write_post_curation_record,
    write_run_manifest,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Run PGA v2.0 on an Oatk plastid FASTA.")
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--annotation-fasta", type=Path, required=True)
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
    return parser.parse_args()


def recreate_dir(path):
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def input_fasta_topology(path):
    topologies = []
    headers = []
    lengths = []
    record_ids = []
    current_length = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current_length is not None:
                    lengths.append(current_length)
                current_length = 0
                header = line[1:].strip()
                record_id = header.split()[0] if header.split() else ""
                if not record_id:
                    raise ValueError(
                        f"Could not determine a FASTA record ID from {path}."
                    )
                headers.append(header)
                record_ids.append(record_id)
                topology = topology_from_fasta_header(header)
                if topology is not None:
                    topologies.append(topology)
            elif current_length is not None:
                current_length += len("".join(line.split()))
    if current_length is not None:
        lengths.append(current_length)

    unique_topologies = sorted(set(topologies))
    topology = unique_topologies[0] if len(unique_topologies) == 1 else None
    sequence_length = lengths[0] if len(lengths) == 1 else None
    record_id = record_ids[0] if len(record_ids) == 1 else None
    return topology, headers, sequence_length, record_id


def run_pga_v2(args, script, reference_dir, target_dir, output_dir, run_root):
    cmd = [
        "perl",
        str(script),
        "-r",
        str(reference_dir),
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
        str(output_dir),
        "-f",
        args.form,
        "-w",
        args.warning,
    ]
    subprocess.run(cmd, check=True, cwd=run_root)
    return cmd


def curate_preliminary_annotation(annotation, record):
    trim_genbank_to_core_sections(annotation)
    curate_genbank_locus(
        annotation,
        topology=record["topology"],
        locus_name=record["id"],
        sequence_length=len(record["sequence"]),
    )
    normalize_genbank_origin_wrapping_locations(annotation)


def main():
    args = parse_args()
    script = args.script.resolve()
    input_fasta = args.input_fasta.resolve()
    annotation_fasta = args.annotation_fasta.resolve()
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
    shutil.copytree(reference_dir, copied_reference_dir)

    input_records = parse_fasta_records(input_fasta)
    circular_origin_rotation = {
        "rotation_applied": False,
        "rotation_offset": 0,
        "new_origin_position": 1,
        "skipped_reason": "input FASTA does not contain exactly one circular record",
    }
    annotation_records = input_records
    preliminary_command = None
    preliminary_annotation = None
    preliminary_selected = None
    if len(input_records) == 1 and input_records[0]["topology"] == "circular":
        preliminary_target_dir = run_root / "preliminary_target"
        preliminary_output_dir = run_root / "preliminary_gb"
        preliminary_target_dir.mkdir(parents=True)
        preliminary_output_dir.mkdir(parents=True)
        preliminary_fasta = preliminary_target_dir / f"{input_fasta.stem}.fa"
        write_fasta_records(input_records, preliminary_fasta)
        preliminary_command = run_pga_v2(
            args,
            script,
            copied_reference_dir,
            preliminary_target_dir,
            preliminary_output_dir,
            run_root,
        )
        preliminary_annotation = run_root / "preliminary.gbk"
        preliminary_selected = copy_first_genbank(
            preliminary_output_dir,
            preliminary_annotation,
        )
        curate_preliminary_annotation(preliminary_annotation, input_records[0])
        circular_origin_rotation = select_circular_origin_rotation_from_genbank(
            preliminary_annotation
        )
        if circular_origin_rotation["rotation_applied"]:
            annotation_records = [
                rotate_fasta_record(
                    input_records[0],
                    circular_origin_rotation["rotation_offset"],
                )
            ]
    elif len(input_records) == 1:
        circular_origin_rotation["skipped_reason"] = (
            f"input record topology is {input_records[0]['topology'] or 'unknown'}"
        )

    write_fasta_records(annotation_records, annotation_fasta)
    target_dir.mkdir(parents=True)
    staged_fasta = target_dir / f"{input_fasta.stem}.fa"
    shutil.copy2(annotation_fasta, staged_fasta)

    cmd = run_pga_v2(
        args,
        script,
        copied_reference_dir,
        target_dir,
        raw_output_dir,
        run_root,
    )
    selected = copy_first_genbank(raw_output_dir, annotation)
    core_sections = trim_genbank_to_core_sections(annotation)
    topology, annotation_headers, sequence_length, input_record_id = input_fasta_topology(
        annotation_fasta
    )
    topology_curation = curate_genbank_locus(
        annotation,
        topology=topology,
        locus_name=input_record_id,
        sequence_length=sequence_length,
    )
    post_curation = curate_genbank_source_metadata(
        annotation,
        args.assembly_name,
        taxid=args.taxid,
        organelle="plastid:chloroplast",
        annotation_method="PGA-Plastid Genome Annotator",
    )
    post_curation["circular_origin_rotation"] = circular_origin_rotation
    post_curation["core_sections"] = core_sections
    post_curation["locus_topology"] = topology_curation
    post_curation["origin_wrapping_locations"] = (
        normalize_genbank_origin_wrapping_locations(annotation)
    )
    post_curation["feature_sort"] = sort_genbank_features_by_location(annotation)
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
            "annotation_fasta": str(annotation_fasta),
            "reference_dir": str(reference_dir),
            "staged_reference_dir": str(copied_reference_dir),
            "input_record_headers": [record["header"] for record in input_records],
            "input_record_id": input_record_id,
            "input_record_length": sequence_length,
            "input_topology": topology,
            "annotation_record_headers": annotation_headers,
            "circular_origin_rotation": circular_origin_rotation,
            "preliminary_command": preliminary_command,
            "preliminary_selected_annotation": (
                str(preliminary_selected) if preliminary_selected else None
            ),
            "preliminary_annotation": (
                str(preliminary_annotation) if preliminary_annotation else None
            ),
            "raw_output_dir": str(raw_output_dir),
            "selected_annotation": str(selected),
            "annotation": str(annotation),
            "post_curation": post_curation,
            **post_curation_record,
            "command": cmd,
        },
    )


if __name__ == "__main__":
    main()
