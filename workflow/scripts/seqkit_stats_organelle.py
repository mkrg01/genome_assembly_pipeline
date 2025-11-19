import subprocess

def run_seqkit_stats(fasta, txt_out, tsv_out, log_err):
    with open(txt_out, "w") as f_txt:
        subprocess.run(
            ["seqkit", "stats", "--all", fasta],
            stdout=f_txt,
            stderr=log_err,
            text=True,
            check=True
        )
    with open(tsv_out, "w") as f_tsv:
        subprocess.run(
            ["seqkit", "stats", "--all", "--tabular", fasta],
            stdout=f_tsv,
            stderr=log_err,
            text=True,
            check=True
        )

def run_seqkit_stats_organelle(snakemake):
    oatk_organelle = snakemake.params.oatk_organelle
    with open(snakemake.log.err, "w") as log_err:
        if oatk_organelle in ("mito", "mito_and_pltd"):
            run_seqkit_stats(
                fasta=str(snakemake.input.mito_ctg_fasta),
                txt_out=str(snakemake.output.mito_txt),
                tsv_out=str(snakemake.output.mito_tsv),
                log_err=log_err
            )
        if oatk_organelle in ("pltd", "mito_and_pltd"):
            run_seqkit_stats(
                fasta=str(snakemake.input.pltd_ctg_fasta),
                txt_out=str(snakemake.output.pltd_txt),
                tsv_out=str(snakemake.output.pltd_tsv),
                log_err=log_err
            )
        if oatk_organelle not in ("mito", "pltd", "mito_and_pltd"):
            raise ValueError(f"Invalid value for 'oatk_organelle' in config.yml: {oatk_organelle}. Must be one of 'mito', 'pltd', or 'mito_and_pltd'.")

run_seqkit_stats_organelle(snakemake)
