import subprocess
from pathlib import Path

def run_cmd(cmd, stdout_file, stderr):
    with open(stdout_file, "w") as fout:
        subprocess.run(cmd, stdout=fout, stderr=stderr, text=True, check=True)

def run_concat_replace(fasta, prefix, out_fasta, stderr):
    cmd_concat = ["seqkit", "concat", fasta, fasta]
    cmd_replace = ["seqkit", "replace", "--pattern", "^", "--replacement", f"{prefix}_"]
    p1 = subprocess.Popen(cmd_concat, stdout=subprocess.PIPE, stderr=stderr, text=True)
    with open(out_fasta, "w") as fout:
        subprocess.run(cmd_replace, stdin=p1.stdout, stdout=fout, stderr=stderr, text=True, check=True)
    p1.stdout.close()
    p1.wait()

def concatenate_organelle_genome(snakemake):
    oatk_organelle = snakemake.params.oatk_organelle
    with open(snakemake.log.err, "w") as log_err:
        if oatk_organelle in ("mito", "mito_and_pltd"):
            run_concat_replace(
                fasta=str(snakemake.input.mito_ctg_fasta),
                prefix="mito",
                out_fasta=str(snakemake.output.mito),
                stderr=log_err,
            )
        if oatk_organelle in ("pltd", "mito_and_pltd"):
            run_concat_replace(
                fasta=str(snakemake.input.pltd_ctg_fasta),
                prefix="pltd",
                out_fasta=str(snakemake.output.pltd),
                stderr=log_err,
            )

        if oatk_organelle == "mito":
            Path(snakemake.output.all_organelle).write_text(
                Path(snakemake.output.mito).read_text()
            )
        elif oatk_organelle == "pltd":
            Path(snakemake.output.all_organelle).write_text(
                Path(snakemake.output.pltd).read_text()
            )
        elif oatk_organelle == "mito_and_pltd":
            with open(snakemake.output.all_organelle, "w") as fout:
                fout.write(Path(snakemake.output.mito).read_text())
                fout.write(Path(snakemake.output.pltd).read_text())
        else:
            raise ValueError(f"Invalid value for 'oatk_organelle' in config.yml: {oatk_organelle}. Must be one of 'mito', 'pltd', or 'mito_and_pltd'.")

concatenate_organelle_genome(snakemake)
