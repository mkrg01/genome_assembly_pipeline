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
    concatenated_fastas = []
    with open(snakemake.log.err, "w") as log_err:
        if "mitochondrion" in oatk_organelle:
            run_concat_replace(
                fasta=str(snakemake.input.mito_ctg_fasta),
                prefix="mitochondrion",
                out_fasta=str(snakemake.output.mito),
                stderr=log_err,
            )
            concatenated_fastas.append(snakemake.output.mito)
        if "chloroplast" in oatk_organelle:
            run_concat_replace(
                fasta=str(snakemake.input.pltd_ctg_fasta),
                prefix="chloroplast",
                out_fasta=str(snakemake.output.pltd),
                stderr=log_err,
            )
            concatenated_fastas.append(snakemake.output.pltd)

        with open(snakemake.output.all_organelle, "w") as fout:
            for fasta in concatenated_fastas:
                fout.write(Path(fasta).read_text())

if __name__ == "__main__":
    concatenate_organelle_genome(snakemake)
