import argparse
import shutil
import subprocess
from pathlib import Path

from organelle_annotation_utils import (
    copy_first_genbank,
    curate_genbank_locus,
    curate_genbank_source_metadata,
    topology_from_fasta_header,
    trim_genbank_to_core_sections,
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
    core_sections = trim_genbank_to_core_sections(annotation)
    topology, input_headers, sequence_length, input_record_id = input_fasta_topology(
        input_fasta
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
    post_curation["core_sections"] = core_sections
    post_curation["locus_topology"] = topology_curation
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
            "input_record_headers": input_headers,
            "input_record_id": input_record_id,
            "input_record_length": sequence_length,
            "input_topology": topology,
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
