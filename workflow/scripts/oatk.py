from pathlib import Path
import subprocess
import sys

def run_oatk(snakemake):
    input = snakemake.input
    output = snakemake.output
    params = snakemake.params
    wildcards = snakemake.wildcards
    threads = snakemake.threads
    outdir = Path(output.utg_final_gfa).parent
    outprefix = outdir / wildcards.assembly_name
    oatk_organelle = params.oatk_organelle

    cmd = ["oatk"]

    if oatk_organelle == "mito":
        cmd += ["-m", str(input.mito_fam)]
    elif oatk_organelle == "pltd":
        cmd += ["-p", str(input.pltd_fam)]
    elif oatk_organelle == "mito_and_pltd":
        cmd += ["-m", str(input.mito_fam), "-p", str(input.pltd_fam)]
    else:
        raise ValueError(f"Invalid value for 'oatk_organelle' in config.yml: {oatk_organelle}. Must be one of 'mito', 'pltd', or 'mito_and_pltd'.")

    cmd += [
        "-o", str(outprefix),
        "-c", str(params.oatk_minimum_kmer_coverage),
        "-t", str(threads),
        str(input.hifi_reads)
    ]

    with open(snakemake.log.out, "w") as log_out, open(snakemake.log.err, "w") as log_err:
        proc = subprocess.run(cmd, stdout=log_out, stderr=log_err)
        if proc.returncode != 0:
            sys.exit(proc.returncode)

run_oatk(snakemake)
