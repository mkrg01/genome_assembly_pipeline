from pathlib import Path
import subprocess
import sys


def read_minimum_kmer_coverage(path):
    value = Path(path).read_text().strip()
    if not value.isdigit() or int(value) <= 0:
        raise ValueError(
            f"Resolved Oatk minimum kmer coverage must be a positive integer, got {value!r}."
        )
    return value


def run_oatk(snakemake):
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    wildcards = snakemake.wildcards
    threads = snakemake.threads
    outdir = Path(output.utg_final_gfa).parent
    outprefix = outdir / wildcards.assembly_name
    oatk_organelle = params.oatk_organelle
    minimum_kmer_coverage = read_minimum_kmer_coverage(input.minimum_kmer_coverage)

    cmd = ["oatk"]

    if oatk_organelle in ("mitochondrion", "mito"):
        cmd += ["-m", str(input.mito_fam)]
    elif oatk_organelle in ("chloroplast", "pltd"):
        cmd += ["-p", str(input.pltd_fam)]
    elif oatk_organelle in ("mitochondrion_and_chloroplast", "mito_and_pltd"):
        cmd += ["-m", str(input.mito_fam), "-p", str(input.pltd_fam)]
    else:
        raise ValueError(
            f"Invalid value for 'oatk_organelle' in config.yml: {oatk_organelle}. "
            "Must be one of 'mitochondrion', 'chloroplast', or "
            "'mitochondrion_and_chloroplast'."
        )

    cmd += [
        "-o", str(outprefix),
        "-c", minimum_kmer_coverage,
        "-t", str(threads),
        str(input.hifi_reads)
    ]

    with open(snakemake.log.out, "w") as log_out, open(snakemake.log.err, "w") as log_err:
        proc = subprocess.run(cmd, stdout=log_out, stderr=log_err)
        if proc.returncode != 0:
            sys.exit(proc.returncode)

run_oatk(snakemake)
