import argparse
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

from organelle_annotation_utils import (
    copy_first_genbank,
    curate_genbank_locus,
    curate_genbank_source_metadata,
    format_cds_auto_translation_validation_errors,
    normalize_genbank_cds_qualifier_order,
    normalize_genbank_feature_location_wrapping,
    normalize_genbank_origin_wrapping_locations,
    normalize_pmga_trans_splicing_qualifiers,
    remove_genbank_feature_qualifiers,
    rotate_fasta_record,
    select_circular_origin_rotation_from_genbank,
    sort_genbank_features_by_location,
    topology_from_fasta_header,
    trim_genbank_to_core_sections,
    validate_genbank_cds_auto_translation,
    write_fasta_records,
    write_post_curation_record,
    write_run_manifest,
)


IMAGE_SUFFIXES = (".sif", ".simg")
FASTA_LINE_WIDTH = 80


def make_fasta_record(header, sequence_lines):
    return {
        "header": header,
        "id": fasta_record_id(header),
        "sequence": "".join(sequence_lines),
        "topology": topology_from_fasta_header(header),
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Run PMGA on an Oatk mitochondrial FASTA.")
    parser.add_argument("--pmga-bundle", type=Path, required=True)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--annotation-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--post-curation", type=Path, required=True)
    parser.add_argument("--db", default="1")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--assembly-name")
    parser.add_argument("--taxid")
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


def parse_fasta_records(path):
    records = []
    header = None
    sequence_lines = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(make_fasta_record(header, sequence_lines))
                header = line[1:].strip()
                if not header:
                    raise ValueError(
                        f"Empty FASTA header found in {path} at line {line_number}."
                    )
                sequence_lines = []
                continue

            if header is None:
                if line.strip():
                    raise ValueError(
                        f"Sequence data found before the first FASTA header in {path} "
                        f"at line {line_number}."
                    )
                continue
            sequence_lines.append("".join(line.split()))

    if header is not None:
        records.append(make_fasta_record(header, sequence_lines))
    if not records:
        raise ValueError(f"No FASTA records found in {path}.")
    return records


def fasta_record_id(header):
    record_id = header.split()[0] if header.split() else ""
    if not record_id:
        raise ValueError(f"Could not determine a FASTA record ID from: {header!r}")
    return record_id


def safe_name(value, default="record"):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._-")
    return safe or default


def unique_pmga_prefixes(prefix, records):
    used = set()
    prefixes = []
    base_prefix = safe_name(prefix, default="pmga")
    for index, record in enumerate(records, start=1):
        record_name = safe_name(record["id"], default=f"record{index:03d}")
        candidate = f"{base_prefix}_{record_name}"
        if candidate in used:
            candidate = f"{candidate}_{index:03d}"
        while candidate in used:
            candidate = f"{base_prefix}_{record_name}_{len(used) + 1:03d}"
        used.add(candidate)
        prefixes.append(candidate)
    return prefixes


def write_single_record_fasta(record, path, record_id=None):
    header = (
        record["header"]
        if record_id is None
        else f"{record_id} {record['header']}"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    sequence = record["sequence"]
    wrapped = "\n".join(textwrap.wrap(sequence, FASTA_LINE_WIDTH))
    if wrapped:
        wrapped = f"{wrapped}\n"
    path.write_text(f">{header}\n{wrapped}")


def run_pmga(runtime, container, cwd, input_fasta, db, prefix, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        runtime,
        "exec",
        "--bind",
        f"{cwd}:{cwd}",
        str(container),
        "mgavas_m",
        "-pid",
        prefix,
        "--in",
        str(input_fasta),
        "--db",
        str(db),
        "--outdir",
        str(output_dir),
    ]
    subprocess.run(cmd, check=True)
    return cmd


def skipped_circular_origin_rotation(record, reason):
    return {
        "rotation_applied": False,
        "rotation_offset": 0,
        "new_origin_position": 1,
        "sequence_length": len(record["sequence"]),
        "skipped_reason": reason,
    }


def curate_preliminary_annotation(annotation, record):
    trim_genbank_to_core_sections(annotation)
    curate_genbank_locus(
        annotation,
        topology=record["topology"],
        locus_name=record["id"],
        sequence_length=len(record["sequence"]),
    )
    normalize_genbank_origin_wrapping_locations(annotation)


def prepare_record_for_annotation(
    *,
    runtime,
    container,
    cwd,
    record,
    db,
    record_prefix,
    preliminary_root,
    record_dir_name,
):
    preliminary_command = None
    preliminary_selected = None
    preliminary_annotation = None
    rotation = skipped_circular_origin_rotation(
        record,
        f"input record topology is {record['topology'] or 'unknown'}",
    )

    if record["topology"] == "circular":
        preliminary_dir = preliminary_root / record_dir_name
        preliminary_input = preliminary_dir / "input.fasta"
        preliminary_output = preliminary_dir / "pmga_raw"
        preliminary_annotation = preliminary_dir / f"{record_prefix}.preliminary.gbk"
        write_single_record_fasta(record, preliminary_input)
        preliminary_command = run_pmga(
            runtime,
            container,
            cwd,
            preliminary_input.resolve(),
            db,
            record_prefix,
            preliminary_output,
        )
        preliminary_selected = copy_first_genbank(
            preliminary_output,
            preliminary_annotation,
        )
        curate_preliminary_annotation(preliminary_annotation, record)
        rotation = select_circular_origin_rotation_from_genbank(preliminary_annotation)

    prepared_record = record
    if rotation["rotation_applied"]:
        prepared_record = rotate_fasta_record(record, rotation["rotation_offset"])

    return {
        "record": prepared_record,
        "circular_origin_rotation": rotation,
        "preliminary_command": preliminary_command,
        "preliminary_selected_annotation": (
            str(preliminary_selected) if preliminary_selected else None
        ),
        "preliminary_annotation": (
            str(preliminary_annotation) if preliminary_annotation else None
        ),
    }


def concatenate_genbank_files(paths, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as out_handle:
        for path in paths:
            text = path.read_text().rstrip()
            if not text:
                raise ValueError(f"Cannot concatenate empty GenBank file: {path}")
            out_handle.write(text)
            out_handle.write("\n")


def main():
    args = parse_args()
    runtime = find_container_runtime()
    container = resolve_pmga_container(args.pmga_bundle)
    input_fasta = args.input_fasta.resolve()
    annotation_fasta = args.annotation_fasta.resolve()
    annotation = args.annotation.resolve()
    manifest = args.manifest.resolve()
    post_curation_path = args.post_curation.resolve()
    raw_output_dir = annotation.parent / "pmga_raw"
    if raw_output_dir.exists():
        shutil.rmtree(raw_output_dir)
    raw_output_dir.mkdir(parents=True, exist_ok=True)
    preliminary_root = annotation.parent / "pmga_preliminary"
    if preliminary_root.exists():
        shutil.rmtree(preliminary_root)
    preliminary_root.mkdir(parents=True, exist_ok=True)

    cwd = Path.cwd().resolve()
    assembly_name = args.assembly_name or args.prefix
    records = parse_fasta_records(input_fasta)
    record_count = len(records)
    record_prefixes = unique_pmga_prefixes(args.prefix, records)
    prepared_records = []
    record_preparations = []
    for index, (record, record_prefix) in enumerate(
        zip(records, record_prefixes),
        start=1,
    ):
        record_dir_name = f"record_{index:03d}_{safe_name(record['id'])}"
        preparation = prepare_record_for_annotation(
            runtime=runtime,
            container=container,
            cwd=cwd,
            record=record,
            db=args.db,
            record_prefix=record_prefix,
            preliminary_root=preliminary_root,
            record_dir_name=record_dir_name,
        )
        prepared_records.append(preparation["record"])
        record_preparations.append(preparation)
    write_fasta_records(prepared_records, annotation_fasta)

    if record_count == 1:
        cmd = run_pmga(
            runtime,
            container,
            cwd,
            annotation_fasta,
            args.db,
            args.prefix,
            raw_output_dir,
        )
        selected = copy_first_genbank(raw_output_dir, annotation)
        core_sections = trim_genbank_to_core_sections(annotation)
        topology_curation = curate_genbank_locus(
            annotation,
            topology=prepared_records[0]["topology"],
            locus_name=prepared_records[0]["id"],
            sequence_length=len(prepared_records[0]["sequence"]),
        )
        post_curation = curate_genbank_source_metadata(
            annotation,
            assembly_name,
            taxid=args.taxid,
            organelle="mitochondrion",
            annotation_method="PMGA-Plant Mitochondrial Genome Annotator",
        )
        post_curation["circular_origin_rotation"] = record_preparations[0][
            "circular_origin_rotation"
        ]
        post_curation["core_sections"] = core_sections
        post_curation["locus_topology"] = topology_curation
        post_curation["feature_location_wrapping"] = (
            normalize_genbank_feature_location_wrapping(annotation)
        )
        post_curation["origin_wrapping_locations"] = (
            normalize_genbank_origin_wrapping_locations(annotation)
        )
        post_curation["pmga_trans_splicing"] = (
            normalize_pmga_trans_splicing_qualifiers(annotation)
        )
        translation_validation = validate_genbank_cds_auto_translation(annotation)
        post_curation["cds_auto_translation_validation"] = translation_validation
        if not translation_validation["cds_auto_translation_passed"]:
            raise RuntimeError(
                format_cds_auto_translation_validation_errors(translation_validation)
            )
        post_curation["feature_qualifier_removal"] = remove_genbank_feature_qualifiers(
            annotation,
            feature_keys=("CDS",),
            qualifiers=("translation", "transl_table", "misc_feature"),
        )
        post_curation["feature_sort"] = sort_genbank_features_by_location(annotation)
        post_curation["cds_qualifier_order"] = normalize_genbank_cds_qualifier_order(
            annotation,
        )
        commands = [cmd]
        selected_annotations = [str(selected)]
        record_manifests = [
            {
                "input_record_id": records[0]["id"],
                "input_record_header": records[0]["header"],
                "input_record_length": len(records[0]["sequence"]),
                "input_record_topology": records[0]["topology"],
                "annotation_record_header": prepared_records[0]["header"],
                "annotation_record_length": len(prepared_records[0]["sequence"]),
                "pmga_prefix": args.prefix,
                "circular_origin_rotation": record_preparations[0][
                    "circular_origin_rotation"
                ],
                "preliminary_command": record_preparations[0][
                    "preliminary_command"
                ],
                "preliminary_selected_annotation": record_preparations[0][
                    "preliminary_selected_annotation"
                ],
                "preliminary_annotation": record_preparations[0][
                    "preliminary_annotation"
                ],
                "raw_output_dir": str(raw_output_dir),
                "selected_annotation": str(selected),
            }
        ]
    else:
        commands = []
        selected_annotations = []
        record_manifests = []
        record_curations = []
        curated_record_annotations = []
        for index, (record, record_prefix) in enumerate(
            zip(prepared_records, record_prefixes),
            start=1,
        ):
            record_dir = raw_output_dir / f"record_{index:03d}_{safe_name(record['id'])}"
            record_input = record_dir / "input.fasta"
            record_pmga_output = record_dir / "pmga_raw"
            record_annotation = record_dir / f"{record_prefix}.gbk"
            original_record = records[index - 1]
            preparation = record_preparations[index - 1]

            write_single_record_fasta(record, record_input)
            cmd = run_pmga(
                runtime,
                container,
                cwd,
                record_input.resolve(),
                args.db,
                record_prefix,
                record_pmga_output,
            )
            selected = copy_first_genbank(record_pmga_output, record_annotation)
            core_sections = trim_genbank_to_core_sections(record_annotation)
            topology_curation = curate_genbank_locus(
                record_annotation,
                topology=record["topology"],
                locus_name=record["id"],
                sequence_length=len(record["sequence"]),
            )
            record_post_curation = curate_genbank_source_metadata(
                record_annotation,
                assembly_name,
                taxid=args.taxid,
                organelle="mitochondrion",
                annotation_method="PMGA-Plant Mitochondrial Genome Annotator",
            )
            record_post_curation["circular_origin_rotation"] = preparation[
                "circular_origin_rotation"
            ]
            record_post_curation["core_sections"] = core_sections
            record_post_curation["locus_topology"] = topology_curation
            record_post_curation["feature_location_wrapping"] = (
                normalize_genbank_feature_location_wrapping(record_annotation)
            )
            record_post_curation["origin_wrapping_locations"] = (
                normalize_genbank_origin_wrapping_locations(record_annotation)
            )
            record_post_curation["pmga_trans_splicing"] = (
                normalize_pmga_trans_splicing_qualifiers(record_annotation)
            )
            translation_validation = validate_genbank_cds_auto_translation(
                record_annotation
            )
            record_post_curation["cds_auto_translation_validation"] = (
                translation_validation
            )
            if not translation_validation["cds_auto_translation_passed"]:
                raise RuntimeError(
                    format_cds_auto_translation_validation_errors(
                        translation_validation
                    )
                )
            record_post_curation["feature_qualifier_removal"] = (
                remove_genbank_feature_qualifiers(
                    record_annotation,
                    feature_keys=("CDS",),
                    qualifiers=("translation", "transl_table", "misc_feature"),
                )
            )
            record_post_curation["feature_sort"] = sort_genbank_features_by_location(
                record_annotation
            )
            record_post_curation["cds_qualifier_order"] = (
                normalize_genbank_cds_qualifier_order(record_annotation)
            )

            commands.append(cmd)
            selected_annotations.append(str(selected))
            curated_record_annotations.append(record_annotation)
            record_manifests.append(
                {
                    "input_record_id": original_record["id"],
                    "input_record_header": original_record["header"],
                    "input_record_length": len(original_record["sequence"]),
                    "input_record_topology": original_record["topology"],
                    "annotation_record_header": record["header"],
                    "annotation_record_length": len(record["sequence"]),
                    "pmga_prefix": record_prefix,
                    "circular_origin_rotation": preparation[
                        "circular_origin_rotation"
                    ],
                    "preliminary_command": preparation["preliminary_command"],
                    "preliminary_selected_annotation": preparation[
                        "preliminary_selected_annotation"
                    ],
                    "preliminary_annotation": preparation["preliminary_annotation"],
                    "split_input_fasta": str(record_input),
                    "raw_output_dir": str(record_pmga_output),
                    "selected_annotation": str(selected),
                    "curated_annotation": str(record_annotation),
                }
            )
            record_curations.append(
                {
                    "input_record_id": original_record["id"],
                    "annotation": str(record_annotation),
                    "post_curation": record_post_curation,
                }
            )

        concatenate_genbank_files(curated_record_annotations, annotation)
        post_curation = {"records": record_curations}

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
            "annotation_fasta": str(annotation_fasta),
            "input_record_count": record_count,
            "records": record_manifests,
            "raw_output_dir": str(raw_output_dir),
            "selected_annotation": selected_annotations[0],
            "selected_annotations": selected_annotations,
            "annotation": str(annotation),
            "post_curation": post_curation,
            **post_curation_record,
            "db": str(args.db),
            "command": commands[0],
            "commands": commands,
        },
    )


if __name__ == "__main__":
    main()
