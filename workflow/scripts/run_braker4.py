"""Run the pinned BRAKER4 workflow on the host and export validated gene models.

The child Snakemake controls Apptainer. Its work directory is deliberately not a
Snakemake directory() output: parent-job retries must not delete its checkpoints.
"""

import argparse
import configparser
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from filter_gff_by_fasta_ids import iter_gff_records, parse_attributes, split_parents


CORE_FILES = ("braker.gtf", "braker.gff3", "braker.codingseq", "braker.aa")
SAMPLE_COLUMNS = (
    "sample_name", "genome", "genome_masked", "protein_fasta", "bam_files",
    "fastq_r1", "fastq_r2", "sra_ids", "varus_genus", "varus_species",
    "isoseq_bam", "isoseq_fastq", "busco_lineage", "reference_gtf",
)


def write_if_changed(path, text):
    path = Path(path)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def checked_path(path):
    path = Path(path).resolve()
    # The upstream shell rules do not consistently quote their path arguments.
    if not re.fullmatch(r"[A-Za-z0-9_./+\-]+", str(path)):
        raise ValueError(f"BRAKER4 requires paths without spaces or shell metacharacters: {path}")
    return path


def file_identity(path):
    path = checked_path(path)
    stat = path.stat()
    if not path.is_file() or stat.st_size == 0:
        raise ValueError(f"BRAKER4 input must be a non-empty file: {path}")
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def prepare_run(args):
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", args.sample):
        raise ValueError("BRAKER4 sample name must start with a letter and contain letters, digits or underscores")
    if not re.fullmatch(r"[A-Za-z0-9_]+_odb12", args.lineage):
        raise ValueError("BRAKER4 v0.5.0-beta requires an odb12 BUSCO lineage, e.g. embryophyta_odb12")
    if args.rnaseq_source == "local":
        if not args.rnaseq_r1 or len(args.rnaseq_r1) != len(args.rnaseq_r2):
            raise ValueError("BRAKER4 local mode requires matching, non-empty R1 and R2 lists")
        if args.varus_genus or args.varus_species:
            raise ValueError("VARUS genus/species cannot be supplied in local RNA-seq mode")
    elif args.rnaseq_source == "varus":
        if args.rnaseq_r1 or args.rnaseq_r2:
            raise ValueError("VARUS mode selects RNA-seq from SRA; do not also supply local FASTQ files")
        if not all(re.fullmatch(r"[A-Za-z][A-Za-z-]*", value or "") for value in (args.varus_genus, args.varus_species)):
            raise ValueError("VARUS requires a valid genus and species")
    else:
        raise ValueError("RNA-seq source must be 'local' or 'varus'")
    if args.threads < 1 or args.mem_mb < 16000:
        raise ValueError("BRAKER4 requires at least 1 thread and 16000 MB memory")
    source = checked_path(args.source)
    if not (source / "Snakefile").is_file():
        raise FileNotFoundError(f"BRAKER4 Snakefile not found in {source}")
    busco = checked_path(args.busco_downloads)
    lineage_dir = busco / "lineages" / args.lineage
    if not (lineage_dir / "dataset.cfg").is_file():
        raise FileNotFoundError(f"Pre-downloaded BUSCO lineage is missing: {lineage_dir}")

    identity = {
        "source": json.loads(Path(args.source_manifest).read_text()),
        "sample": args.sample, "lineage": args.lineage,
        "rnaseq_source": args.rnaseq_source,
        "varus_genus": args.varus_genus, "varus_species": args.varus_species,
        "genome": file_identity(args.genome), "proteins": file_identity(args.proteins),
        "rnaseq_r1": [file_identity(p) for p in args.rnaseq_r1],
        "rnaseq_r2": [file_identity(p) for p in args.rnaseq_r2],
        "template_sha256": hashlib.sha256(Path(args.template).read_bytes()).hexdigest(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    # A changed input set or model configuration gets a fresh run, while retries
    # with the same inputs reuse all child metadata. Resource changes can resume.
    run_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
    run_dir = checked_path(args.output_dir) / "work" / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(exist_ok=True)
    r1, r2 = [], []
    for number, pair in enumerate(zip(identity["rnaseq_r1"], identity["rnaseq_r2"]), 1):
        for mate, item, paths in zip((1, 2), pair, (r1, r2)):
            # Upstream truncates IDs at the first dot. Unique aliases preserve
            # distinct libraries such as tissue.rep1 and tissue.rep2.
            link = inputs_dir / f"rnaseq_{number:04d}_{mate}.fastq"
            if not link.is_symlink():
                link.symlink_to(item["path"])
            paths.append(str(link))
    row = dict.fromkeys(SAMPLE_COLUMNS, "")
    row.update(sample_name=args.sample, genome=identity["genome"]["path"],
               genome_masked=identity["genome"]["path"], protein_fasta=identity["proteins"]["path"],
               fastq_r1=":".join(r1), fastq_r2=":".join(r2), busco_lineage=args.lineage,
               varus_genus=args.varus_genus, varus_species=args.varus_species)
    csv_text = io.StringIO()
    writer = csv.DictWriter(csv_text, fieldnames=SAMPLE_COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerow(row)
    write_if_changed(run_dir / "samples.csv", csv_text.getvalue())

    ini = configparser.ConfigParser()
    ini.read(args.template)
    ini.set("PARAMS", "use_varus", "1" if args.rnaseq_source == "varus" else "0")
    ini["paths"] = {
        "samples_file": str(run_dir / "samples.csv"),
        "augustus_config_path": str(run_dir / "augustus_config"),
        "busco_download_path": str(busco),
        "compleasm_download_path": str(busco / "lineages"),
    }
    ini["SLURM_ARGS"] = {
        "cpus_per_task": str(args.threads), "mem_of_node": str(args.mem_mb),
        "max_runtime": "4320",
    }
    ini_text = io.StringIO()
    ini.write(ini_text)
    write_if_changed(run_dir / "config.ini", ini_text.getvalue())
    write_if_changed(run_dir / "inputs.json", json.dumps(identity, indent=2) + "\n")

    cache = checked_path(args.container_cache)
    cache.mkdir(parents=True, exist_ok=True)
    tmpdir = run_dir / "tmp"
    tmpdir.mkdir(exist_ok=True)
    bind_paths = {str(run_dir), str(source), str(busco), str(cache)}
    for item in [identity["genome"], identity["proteins"], *identity["rnaseq_r1"], *identity["rnaseq_r2"]]:
        bind_paths.add(str(Path(item["path"]).parent))
    command = [
        sys.executable, "-m", "snakemake", "--snakefile", str(source / "Snakefile"),
        "--directory", str(run_dir), "--executor", "local", "--scheduler", "greedy",
        "--profile", "none", "--workflow-profile", "none",
        "--cores", str(args.threads), "--resources", f"mem_mb={args.mem_mb}",
        "--default-resources", f"tmpdir={tmpdir}",
        "--use-singularity", "--singularity-prefix", str(cache),
        "--singularity-args", shlex.join(["--bind", ",".join(sorted(bind_paths))]),
        "--rerun-incomplete", "--latency-wait", "60", "--printshellcmds", "all",
    ]
    if args.dry_run:
        command.append("--dry-run")
    collected = run_dir / "output" / args.sample / "results"
    if (collected / ".done").exists() and any(not (collected / (name + ".gz")).is_file() for name in CORE_FILES):
        # Upstream declares only .done, not the compressed files it collects.
        command.extend(["--forcerun", "collect_results"])
    # The generated config is authoritative. Inherited BRAKER4_* variables and
    # cluster profiles must not change the model or submit nested scheduler jobs.
    env = {key: value for key, value in os.environ.items() if not key.startswith("BRAKER4_")}
    env.update(BRAKER4_CONFIG=str(run_dir / "config.ini"), TMPDIR=str(tmpdir),
               APPTAINERENV_TMPDIR=str(tmpdir), SINGULARITYENV_TMPDIR=str(tmpdir))
    return run_dir, identity, command, env


def fasta_ids(path):
    ids = set()
    current = None
    has_sequence = False
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current is not None and not has_sequence:
                    raise ValueError(f"Empty sequence for {current} in {path}")
                fields = line[1:].split()
                if not fields or fields[0] in ids:
                    raise ValueError(f"Missing or duplicate FASTA ID in {path}: {line.strip()}")
                current = fields[0]
                ids.add(current)
                has_sequence = False
            elif line.strip():
                if current is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                has_sequence = True
    if not ids or not has_sequence:
        raise ValueError(f"Empty FASTA or final sequence in {path}")
    return ids


def validate_outputs(directory):
    cds_ids = fasta_ids(directory / "braker.codingseq")
    if cds_ids != fasta_ids(directory / "braker.aa"):
        raise ValueError("BRAKER4 CDS and protein FASTA IDs differ")
    genes, transcripts, cds_parents = set(), {}, set()
    for _, _, fields in iter_gff_records(directory / "braker.gff3"):
        if fields is None:
            continue
        attrs = parse_attributes(fields[8])
        if fields[2] == "gene":
            gene = attrs.get("ID")
            if not gene or gene in genes:
                raise ValueError("Missing or duplicate gene ID in BRAKER4 GFF3")
            genes.add(gene)
        elif fields[2] == "mRNA":
            tx = attrs.get("ID")
            parents = split_parents(attrs.get("Parent", ""))
            if not tx or tx in transcripts or len(parents) != 1:
                raise ValueError("Invalid mRNA ID/Parent in BRAKER4 GFF3")
            transcripts[tx] = parents[0]
        elif fields[2] == "CDS":
            cds_parents.update(split_parents(attrs.get("Parent", "")))
    if set(transcripts) != cds_ids or cds_parents != cds_ids:
        raise ValueError("BRAKER4 GFF3 mRNA/CDS IDs do not match the FASTA IDs")
    if set(transcripts.values()) != genes:
        raise ValueError("BRAKER4 GFF3 gene/mRNA relationships are incomplete")
    gtf_ids = set()
    with (directory / "braker.gtf").open() as handle:
        for line in handle:
            if not line.startswith("#"):
                gtf_ids.update(re.findall(r'transcript_id "([^"]+)"', line))
    if gtf_ids != cds_ids:
        raise ValueError("BRAKER4 GTF transcript IDs do not match the FASTA IDs")
    return {"genes": len(genes), "transcripts": len(transcripts)}


def export_results(run_dir, output_dir, sample, identity, command):
    results = run_dir / "output" / sample / "results"
    if not (results / ".done").is_file():
        raise FileNotFoundError(f"BRAKER4 did not finish results collection: {results}")
    output_dir = Path(output_dir)
    with tempfile.TemporaryDirectory(dir=output_dir) as tmp:
        staged = Path(tmp)
        for name in CORE_FILES:
            with gzip.open(results / (name + ".gz"), "rb") as source, (staged / name).open("wb") as out:
                shutil.copyfileobj(source, out)
        counts = validate_outputs(staged)
        manifest = {**identity, "workdir": str(run_dir), "results": str(results),
                    "command": command, "counts": counts}
        (staged / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
        for name in (*CORE_FILES, "run.json"):
            (staged / name).replace(output_dir / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "source-manifest", "template", "genome", "proteins",
                 "busco-downloads", "container-cache", "output-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--rnaseq-source", choices=("local", "varus"), default="local")
    parser.add_argument("--varus-genus", default="")
    parser.add_argument("--varus-species", default="")
    parser.add_argument("--rnaseq-r1", nargs="*", type=Path, default=[])
    parser.add_argument("--rnaseq-r2", nargs="*", type=Path, default=[])
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--mem-mb", type=int, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    run_dir, identity, command, env = prepare_run(args)
    print(f"BRAKER4 work directory: {run_dir}", flush=True)
    print(shlex.join(command), flush=True)
    subprocess.run(command, env=env, check=True)
    if not args.dry_run:
        export_results(run_dir, args.output_dir, args.sample, identity, command)


if __name__ == "__main__":
    main()
