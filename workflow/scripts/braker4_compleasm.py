"""Prepare a run-private Compleasm cache and verify BRAKER4 rescue outputs."""

import configparser
import fcntl
import gzip
import hashlib
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.request import urlopen


BUSCO_CATALOG_URL = "https://busco-data.ezlab.org/v5/data/file_versions.tsv"


def complete_summary(path):
    text = Path(path).read_text()
    for key in ("S", "D", "F", "M"):
        if not re.search(rf"(?m)^\s*{key}:\s*\d", text):
            raise ValueError(f"Incomplete Compleasm summary: {path}")


def private_cache(run_dir):
    return Path(run_dir) / "compleasm_downloads"


def container_command(args, run_dir):
    ini = configparser.ConfigParser()
    ini.read(args.template)
    uri = ini["containers"]["braker3_image"]
    cache = Path(args.container_cache).resolve()
    # Match the name used by the child Snakemake workflow.
    image = cache / (hashlib.md5(uri.encode()).hexdigest() + ".simg")
    if not image.is_file():
        cache.mkdir(parents=True, exist_ok=True)
        with image.with_suffix(".pull.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if not image.exists():
                with tempfile.TemporaryDirectory(dir=cache, prefix=".braker3-pull-") as tmp:
                    staged = Path(tmp) / image.name
                    subprocess.run(["apptainer", "pull", str(staged), uri], check=True, timeout=3600)
                    if not staged.is_file() or not staged.stat().st_size:
                        raise ValueError(f"Apptainer did not create the BRAKER3 image for {uri}")
                    staged.rename(image)
    if not image.is_file() or not image.stat().st_size:
        raise ValueError(f"Invalid BRAKER3 image in cache: {image}")
    return ["apptainer", "exec", "--cleanenv", "--no-home", "--bind", str(run_dir),
            str(image), "python3"]


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verify_lineage(source, candidate):
    profiles = {p.name for p in (source / "hmms").glob("*.hmm")}
    if not profiles or profiles != {p.name for p in (candidate / "hmms").iterdir()}:
        raise ValueError(f"Compleasm and BUSCO HMM profile inventories differ: {candidate}")
    paths = ["scores_cutoff", "refseq_db.faa.gz", *(f"hmms/{name}" for name in sorted(profiles))]
    # The pinned Compleasm downloader only extracts the reference proteins,
    # HMMs and score cutoffs. BRAKER4 accepts odb12 datasets without lengths.
    if (candidate / "lengths_cutoff").exists():
        paths.append("lengths_cutoff")
    for name in paths:
        if digest(source / name) != digest(candidate / name):
            raise ValueError(f"Compleasm lineage differs from the selected BUSCO dataset: {name}")


def write_fixed_lineage_catalog(cache, lineage, catalog, checksum):
    if hashlib.md5(catalog).hexdigest() != checksum:
        raise ValueError("BUSCO catalog checksum does not match")
    required = {lineage, "eukaryota_odb12"}
    rows = {}
    for line in catalog.decode().splitlines():
        name, date, md5, _, kind = line.split()
        if kind == "lineages" and name in required:
            if name in rows or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date) or not re.fullmatch(r"[a-f0-9]{32}", md5):
                raise ValueError(f"Invalid BUSCO catalog entry for {name}")
            rows[name] = line
    if rows.keys() != required:
        raise ValueError("Fixed lineages are missing from the BUSCO catalog")

    # BRAKER4 names an explicit lineage in all Compleasm calls. The pinned
    # downloader cannot parse newer auto-lineage placement IDs (e.g. odb12.2).
    # Install a private catalog with only the fixed lineage datasets it uses.
    data = ("\n".join(rows.values()) + "\n").encode()
    cache.mkdir()
    (cache / "file_versions.upstream.tsv").write_bytes(catalog)
    (cache / "file_versions.upstream.tsv.hash").write_text(checksum + "\n")
    (cache / "file_versions.tsv").write_bytes(data)
    (cache / "file_versions.tsv.hash").write_text(hashlib.md5(data).hexdigest() + "  file_versions.tsv\n")
    (cache / "file_versions.tsv.done").touch()


def download_cache(args, run_dir, command):
    cache = private_cache(run_dir)
    source = Path(args.busco_downloads) / "lineages" / args.lineage
    if cache.exists():
        if list(cache.glob("*.tmp")):
            raise ValueError(f"Incomplete Compleasm download requires review: {cache}")
        verify_lineage(source, cache / args.lineage)
        return cache
    with tempfile.TemporaryDirectory(dir=run_dir, prefix=".compleasm-cache-") as tmp:
        staged = Path(tmp) / "library"
        with urlopen(BUSCO_CATALOG_URL + ".hash", timeout=60) as response:
            checksum = response.read().decode().split()[0]
        with urlopen(BUSCO_CATALOG_URL, timeout=60) as response:
            catalog = response.read()
        write_fixed_lineage_catalog(staged, args.lineage, catalog, checksum)
        with (Path(run_dir) / "compleasm_download.log").open("w") as log:
            subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "download",
                            args.lineage, "-L", str(staged)], check=True, stdout=log,
                           stderr=subprocess.STDOUT, timeout=3600)
        lineage = staged / args.lineage
        verify_lineage(source, lineage)
        if (source / "lengths_cutoff").is_file() and not (lineage / "lengths_cutoff").exists():
            shutil.copyfile(source / "lengths_cutoff", lineage / "lengths_cutoff")
        staged.rename(cache)
    return cache


def write_probe(reference, probe):
    with gzip.open(reference, "rt") as source, probe.open("w") as output:
        header = next(source, "")
        if not header.startswith(">"):
            raise ValueError(f"Compleasm reference is not a FASTA: {reference}")
        output.write(header)
        for line in source:
            if line.startswith(">"):
                break
            output.write(line)


def prepare_cache(args, run_dir):
    """Use the pinned container's unmodified downloader before admitting a run."""
    command = container_command(args, run_dir)
    cache = download_cache(args, run_dir, command)
    # Exercise the same protein-mode startup used by best_by_compleasm.
    reference = cache / args.lineage / "refseq_db.faa.gz"
    probe = Path(run_dir) / "compleasm_preflight.faa"
    write_probe(reference, probe)
    with tempfile.TemporaryDirectory(dir=run_dir, prefix=".compleasm-probe-") as tmp:
        output = Path(tmp) / "results"
        with (Path(run_dir) / "compleasm_preflight.log").open("w") as log:
            subprocess.run([*command, "/opt/compleasm_kit/compleasm.py", "protein",
                            "-p", str(probe), "-l", args.lineage, "-L", str(cache),
                            "-o", str(output), "-t", "1"],
                           check=True, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        complete_summary(output / "summary.txt")


def validate_rescue(run_dir, sample):
    root = Path(run_dir) / "output" / sample
    log = (root / "best_by_compleasm.log").read_text()
    if re.search(r"Return code of subprocess was [1-9]|Traceback|script may have failed", log):
        raise ValueError("BRAKER4 best_by_compleasm failed; exported predictions are not accepted")
    # Upstream suppresses subprocess failures. A result for every compared set
    # is required even when the original TSEBRA prediction remains the best.
    for name in ("genemark", "augustus", "braker"):
        folder = root / "best_by_compleasm_tmp" / "single" / name
        summaries = list(folder.rglob("summary.txt"))
        if len(summaries) != 1:
            raise ValueError(f"Missing Compleasm rescue evaluation for {name}: {folder}")
        complete_summary(summaries[0])
    complete_summary(root / "compleasm_proteins" / "summary.txt")
