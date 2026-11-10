circos_ids = [track["id"] for track in config["circos_plot_tracks"]]
circos_id2tracks = {track["id"]: track for track in config["circos_plot_tracks"]}
repeat_classes = [circos_id for circos_id in circos_ids if circos_id in ["LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"]]

wildcard_constraints:
    assembly_name = organism_name,
    selected_assembly = selected_assembly_pattern,
    repeat_class = "|".join(["LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"])

rule calculate_long_contig_lengths:
    input:
        f"results/{downstream_assembly_name}/assembly_long_contigs/{{selected_assembly}}/{{assembly_name}}.fa"
    output:
        "results/circos_plot/{selected_assembly}/{assembly_name}_long_contig_length.tsv"
    log:
        out = "logs/calculate_long_contig_lengths_{selected_assembly}_{assembly_name}.out",
        err = "logs/calculate_long_contig_lengths_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        (
            seqkit fx2tab --name --length {input:q} > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule make_gene_bed:
    input:
        lambda w: external_gff3 if external_gff3 is not None else (
            f"results/braker4/{w.selected_assembly}/{w.assembly_name}/braker.gff3"
        )
    output:
        "results/circos_plot/gene/{selected_assembly}/{assembly_name}_gene.bed"
    log:
        out = "logs/make_gene_bed_{selected_assembly}_{assembly_name}.out",
        err = "logs/make_gene_bed_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        (
            grep -v '^#' {input:q} \
            | awk '$3=="gene"{{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}}' \
            > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule calculate_gene_coverage_per_window:
    input:
        gene_bed = "results/circos_plot/gene/{selected_assembly}/{assembly_name}_gene.bed",
        long_contigs = "results/circos_plot/{selected_assembly}/{assembly_name}_long_contig_length.tsv"
    output:
        windows = "results/circos_plot/gene/{selected_assembly}/{assembly_name}_windows.bed",
        windows_gene_coverage = "results/circos_plot/gene/{selected_assembly}/{assembly_name}_windows_gene_coverage.bed"
    log:
        out = "logs/calculate_gene_coverage_per_window_{selected_assembly}_{assembly_name}.out",
        err = "logs/calculate_gene_coverage_per_window_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    params:
        window_size = lambda w: circos_id2tracks["gene"]["window_size"]
    shell:
        """
        (
            bedtools makewindows \
                -g {input.long_contigs:q} \
                -w {params.window_size} \
                > {output.windows:q}

            bedtools coverage \
                -a {output.windows:q} \
                -b {input.gene_bed:q} \
                > {output.windows_gene_coverage:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule make_repeat_bed:
    input:
        lambda w: external_repeatmasker if external_repeatmasker is not None else (
            f"results/repeatmasker/{w.selected_assembly}/{w.assembly_name}.fa.out.xm"
        )
    output:
        "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{assembly_name}_repeat_all.bed"
    log:
        out = "logs/make_repeat_bed_{selected_assembly}_{assembly_name}.out",
        err = "logs/make_repeat_bed_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    params:
        script = os.path.join(workflow.basedir, "scripts/make_repeat_bed.awk")
    shell:
        """
        (
            awk -f {params.script:q} {input:q} > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule make_repeat_bed_per_class:
    input:
        "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{assembly_name}_repeat_all.bed"
    output:
        "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/{assembly_name}_repeat.bed"
    log:
        out = "logs/make_repeat_bed_per_class_{selected_assembly}_{assembly_name}_{repeat_class}.out",
        err = "logs/make_repeat_bed_per_class_{selected_assembly}_{assembly_name}_{repeat_class}.err"
    conda:
        "../envs/pybase.yml"
    script:
        "../scripts/extract_repeat_class.py"

rule calculate_repeat_coverage_per_window:
    input:
        repeat_bed = "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/{assembly_name}_repeat.bed",
        long_contigs = "results/circos_plot/{selected_assembly}/{assembly_name}_long_contig_length.tsv"
    output:
        windows = "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/{assembly_name}_windows.bed",
        windows_repeat_coverage = "results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/{assembly_name}_windows_repeat_coverage.bed"
    log:
        out = "logs/calculate_repeat_coverage_per_window_{selected_assembly}_{assembly_name}_{repeat_class}.out",
        err = "logs/calculate_repeat_coverage_per_window_{selected_assembly}_{assembly_name}_{repeat_class}.err"
    conda:
        "../envs/bedtools.yml"
    params:
        window_size = lambda wildcards: circos_id2tracks[wildcards.repeat_class]["window_size"]
    shell:
        """
        (
            bedtools makewindows \
                -g {input.long_contigs:q} \
                -w {params.window_size} \
                > {output.windows:q}

            bedtools coverage \
                -a {output.windows:q} \
                -b {input.repeat_bed:q} \
                > {output.windows_repeat_coverage:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule tidk_search_for_circos_plot:
    input:
        f"results/{downstream_assembly_name}/assembly_long_contigs/{{selected_assembly}}/{{assembly_name}}.fa"
    output:
        "results/circos_plot/tidk_repeat/{selected_assembly}/{assembly_name}_windows_tidk.tsv"
    log:
        out = "logs/tidk_search_for_circos_plot_{selected_assembly}_{assembly_name}.out",
        err = "logs/tidk_search_for_circos_plot_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/tidk.yml"
    params:
        telomeric_repeat_unit = config["tidk_telomeric_repeat_unit"],
        window_size = lambda w: circos_id2tracks["tidk"]["window_size"]
    shell:
        """
        (
            tidk search \
                --string {params.telomeric_repeat_unit} \
                --window {params.window_size} \
                --dir $(dirname {output:q}) \
                --output {wildcards.assembly_name} \
                {input:q}
            mv $(dirname {output:q})/{wildcards.assembly_name}_telomeric_repeat_windows.tsv {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """
rule circos_plot:
    input:
        unpack(circos_plot_input_path)
    output:
        "results/circos_plot/{selected_assembly}/{assembly_name}_circos_plot.pdf"
    log:
        out = "logs/circos_plot_{selected_assembly}_{assembly_name}.out",
        err = "logs/circos_plot_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/pycirclize.yml"
    script:
        "../scripts/circos_plot.py"

rule linear_plot:
    input:
        unpack(circos_plot_input_path)
    output:
        "results/circos_plot/{selected_assembly}/{assembly_name}_linear_plot.pdf"
    log:
        out = "logs/linear_plot_{selected_assembly}_{assembly_name}.out",
        err = "logs/linear_plot_{selected_assembly}_{assembly_name}.err"
    conda:
        "../envs/pycirclize.yml"
    script:
        "../scripts/linear_plot.py"
