rule minimap2_self_dotplot:
    input:
        "results/{assembly}/assembly_long_contigs/{selected_assembly}/{assembly_name}.fa"
    output:
        "results/{assembly}/dotplot/{selected_assembly}/{assembly_name}.self.paf"
    log:
        out = "logs/minimap2_self_dotplot_{assembly}_{selected_assembly}_{assembly_name}.out",
        err = "logs/minimap2_self_dotplot_{assembly}_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/minimap2.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        """
        (
            minimap2 -t {threads} -DP -k19 -w19 -m200 \
                {input:q} {input:q} > {output:q}
        ) > {log.out} 2> {log.err}
        """

rule plot_self_dotplot:
    input:
        assembly = "results/{assembly}/assembly_long_contigs/{selected_assembly}/{assembly_name}.fa",
        paf = "results/{assembly}/dotplot/{selected_assembly}/{assembly_name}.self.paf",
        plot_script = os.path.join(workflow.basedir, "scripts", "plot_self_dotplot.py")
    output:
        pdf = f"{SELF_DOTPLOT_STEM}.pdf",
        contigs = f"{SELF_DOTPLOT_STEM}.contigs.tsv",
        summary = f"{SELF_DOTPLOT_STEM}.summary.json",
        additional_pdfs = self_dotplot_additional_paths(SELF_DOTPLOT_STEM, ".pdf"),
        additional_summaries = self_dotplot_additional_paths(SELF_DOTPLOT_STEM, ".summary.json"),
        comparison = self_dotplot_comparison_paths(SELF_DOTPLOT_STEM)
    log:
        out = "logs/plot_self_dotplot_{assembly}_{selected_assembly}_{assembly_name}.out",
        err = "logs/plot_self_dotplot_{assembly}_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/pybase.yml"
    params:
        min_long_contig_length = config["min_long_contig_length"],
        min_alignment_length = self_dotplot_min_alignment_length,
        additional_plots = lambda wildcards, output: [
            argument
            for length, pdf, summary in zip(
                self_dotplot_additional_lengths, output.additional_pdfs, output.additional_summaries,
            )
            for argument in ("--additional-plot", str(length), pdf, summary)
        ],
        comparison_args = lambda wildcards, output: [
            argument for pdf in output.comparison for argument in ("--comparison-pdf", pdf)
        ]
    shell:
        """
        (
            python3 {input.plot_script:q} \
                --fasta {input.assembly:q} \
                --paf {input.paf:q} \
                --pdf {output.pdf:q} \
                --contigs {output.contigs:q} \
                --summary {output.summary:q} \
                --organism-name {wildcards.assembly_name} \
                --min-long-contig-length {params.min_long_contig_length} \
                --min-alignment-length {params.min_alignment_length} \
                {params.additional_plots:q} {params.comparison_args:q}
        ) > {log.out} 2> {log.err}
        """
