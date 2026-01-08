circos_ids = [track["id"] for track in config["circos_plot"]]
circos_id2tracks = {track["id"]: track for track in config["circos_plot"]}
repeat_classes = [circos_id for circos_id in circos_ids if circos_id in ["LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"]]

wildcard_constraints:
    assembly_name = config["assembly_name"],
    repeat_class = "|".join(["LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"])

rule calculate_long_contig_lengths:
    input:
        "results/fcs/assembly_long_contigs/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        "results/circos_plot/{assembly_name}_long_contig_length.tsv"
    log:
        out = "logs/calculate_long_contig_lengths_{assembly_name}.out",
        err = "logs/calculate_long_contig_lengths_{assembly_name}.err"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        (
            seqkit fx2tab --name --length {input} > {output}
        ) > {log.out} 2> {log.err}
        """

rule make_gene_bed:
    input:
        "results/braker3/{assembly_name}/braker.gff3"
    output:
        "results/circos_plot/gene/{assembly_name}_gene.bed"
    log:
        out = "logs/make_gene_bed_{assembly_name}.out",
        err = "logs/make_gene_bed_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        (
            grep -v '^#' {input} \
            | awk '$3=="gene"{{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}}' \
            > {output}
        ) > {log.out} 2> {log.err}
        """

rule calculate_gene_coverage_per_window:
    input:
        gene_bed = "results/circos_plot/gene/{assembly_name}_gene.bed",
        long_contigs = "results/circos_plot/{assembly_name}_long_contig_length.tsv"
    output:
        windows = "results/circos_plot/gene/{assembly_name}_windows.bed",
        windows_gene_coverage = "results/circos_plot/gene/{assembly_name}_windows_gene_coverage.bed"
    log:
        out = "logs/calculate_gene_coverage_per_window_{assembly_name}.out",
        err = "logs/calculate_gene_coverage_per_window_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    params:
        window_size = circos_id2tracks["gene"]["window_size"]
    shell:
        """
        (
            bedtools makewindows \
                -g {input.long_contigs} \
                -w {params.window_size} \
                > {output.windows}

            bedtools coverage \
                -a {output.windows} \
                -b {input.gene_bed} \
                > {output.windows_gene_coverage}
        ) > {log.out} 2> {log.err}
        """

rule make_repeat_bed:
    input:
        "results/repeatmasker/{assembly_name}.asm.bp.p_ctg.fa.out.xm"
    output:
        "results/circos_plot/repeatmasker_repeat/{assembly_name}_repeat_all.bed"
    log:
        out = "logs/make_repeat_bed_{assembly_name}.out",
        err = "logs/make_repeat_bed_{assembly_name}.err"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        (
            awk -f workflow/scripts/make_repeat_bed.awk {input} > {output}
        ) > {log.out} 2> {log.err}
        """

rule make_repeat_bed_per_class:
    input:
        "results/circos_plot/repeatmasker_repeat/{assembly_name}_repeat_all.bed"
    output:
        "results/circos_plot/repeatmasker_repeat/{repeat_class}/{assembly_name}_repeat.bed"
    log:
        out = "logs/make_repeat_bed_per_class_{assembly_name}_{repeat_class}.out",
        err = "logs/make_repeat_bed_per_class_{assembly_name}_{repeat_class}.err"
    conda:
        "../envs/pybase.yml"
    script:
        "../scripts/extract_repeat_class.py"

rule calculate_repeat_coverage_per_window:
    input:
        repeat_bed = "results/circos_plot/repeatmasker_repeat/{repeat_class}/{assembly_name}_repeat.bed",
        long_contigs = "results/circos_plot/{assembly_name}_long_contig_length.tsv"
    output:
        windows = "results/circos_plot/repeatmasker_repeat/{repeat_class}/{assembly_name}_windows.bed",
        windows_repeat_coverage = "results/circos_plot/repeatmasker_repeat/{repeat_class}/{assembly_name}_windows_repeat_coverage.bed"
    log:
        out = "logs/calculate_repeat_coverage_per_window_{assembly_name}_{repeat_class}.out",
        err = "logs/calculate_repeat_coverage_per_window_{assembly_name}_{repeat_class}.err"
    conda:
        "../envs/bedtools.yml"
    params:
        window_size = lambda wildcards: circos_id2tracks[wildcards.repeat_class]["window_size"]
    shell:
        """
        (
            bedtools makewindows \
                -g {input.long_contigs} \
                -w {params.window_size} \
                > {output.windows}

            bedtools coverage \
                -a {output.windows} \
                -b {input.repeat_bed} \
                > {output.windows_repeat_coverage}
        ) > {log.out} 2> {log.err}
        """

rule tidk_search_for_circos_plot:
    input:
        "results/fcs/assembly_long_contigs/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        "results/circos_plot/tidk_repeat/{assembly_name}_windows_tidk.tsv"
    log:
        out = "logs/tidk_search_for_circos_plot_{assembly_name}.out",
        err = "logs/tidk_search_for_circos_plot_{assembly_name}.err"
    conda:
        "../envs/tidk.yml"
    params:
        telomeric_repeat_unit = config["tidk_telomeric_repeat_unit"],
        window_size = circos_id2tracks["tidk"]["window_size"]
    shell:
        """
        (
            tidk search \
                --string {params.telomeric_repeat_unit} \
                --window {params.window_size} \
                --dir $(dirname {output}) \
                --output {wildcards.assembly_name} \
                {input}
            mv $(dirname {output})/{wildcards.assembly_name}_telomeric_repeat_windows.tsv {output}
        ) > {log.out} 2> {log.err}
        """

def circos_plot_input_path():
    input_dict = {"contig": f"results/circos_plot/{{assembly_name}}_long_contig_length.tsv"}
    if "gene" in circos_ids:
        input_dict["gene"] = f"results/circos_plot/gene/{{assembly_name}}_windows_gene_coverage.bed"
    for repeat_class in repeat_classes:
        input_dict[f"{repeat_class}"] = f"results/circos_plot/repeatmasker_repeat/{repeat_class}/{{assembly_name}}_windows_repeat_coverage.bed"
    if "tidk" in circos_ids:
        input_dict["tidk"] = f"results/circos_plot/tidk_repeat/{{assembly_name}}_windows_tidk.tsv"
    return input_dict

rule circos_plot:
    input:
        **circos_plot_input_path()
    output:
        "results/circos_plot/{assembly_name}_circos_plot.pdf"
    log:
        out = "logs/circos_plot_{assembly_name}.out",
        err = "logs/circos_plot_{assembly_name}.err"
    conda:
        "../envs/pycirclize.yml"
    script:
        "../scripts/circos_plot.py"
