rule bam2fastq:
    input:
        "raw_data/{sample_id}.hifi_reads.bam"
    output:
        "results/raw_reads/{sample_id}_hifi_reads.fastq.gz"
    log:
        out = "logs/bam2fastq_{sample_id}.out",
        err = "logs/bam2fastq_{sample_id}.err"
    conda:
        "../envs/bam2fastq.yml"
    threads:
        max(1, int(workflow.cores * 0.5))
    shell:
        "bam2fastq \
            {input} \
            --output $(dirname {output})/$(basename {output} .fastq.gz) \
            --num-threads {threads} > {log.out} 2> {log.err}"

rule fastplong:
    input:
        "results/raw_reads/{sample_id}_hifi_reads.fastq.gz"
    output:
        reads = "results/fastplong/{sample_id}_hifi_reads_curated.fastq.gz",
        report_html = "results/fastplong/{sample_id}_report.html",
        report_json = "results/fastplong/{sample_id}_report.json"
    log:
        out = "logs/fastplong_{sample_id}.out",
        err = "logs/fastplong_{sample_id}.err"
    conda:
        "../envs/fastplong.yml"
    threads:
        max(1, int(workflow.cores * 0.5))
    shell:
        "fastplong \
            --in {input} \
            --out {output.reads} \
            --html {output.report_html} \
            --json {output.report_json} \
            --thread {threads} > {log.out} 2> {log.err}"

rule jellyfish_count:
    input:
        "results/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        "results/jellyfish/{sample_id}_mer_counts.jf"
    log:
        out = "logs/jellyfish_count_{sample_id}.out",
        err = "logs/jellyfish_count_{sample_id}.err"
    conda:
        "../envs/jellyfish.yml"
    threads:
        1
    shell:
        "jellyfish count \
            --mer-len 21 \
            --size 1G \
            --threads {threads} \
            --canonical \
            --output {output} \
            <(zcat {input}) > {log.out} 2> {log.err}"

rule jellyfish_histo:
    input:
        "results/jellyfish/{sample_id}_mer_counts.jf"
    output:
        "results/jellyfish/{sample_id}_jellyfish.histo"
    log:
        out = "logs/jellyfish_histo_{sample_id}.out",
        err = "logs/jellyfish_histo_{sample_id}.err"
    conda:
        "../envs/jellyfish.yml"
    threads:
        1
    shell:
        "jellyfish histo \
            {input} \
            --output {output} \
            --threads {threads} > {log.out} 2> {log.err}"

rule genomescope2:
    input:
        "results/jellyfish/{sample_id}_jellyfish.histo"
    output:
        linear_plot = "results/genomescope2/{sample_id}/linear_plot.png",
        log_plot = "results/genomescope2/{sample_id}/log_plot.png",
        model = "results/genomescope2/{sample_id}/model.txt",
        progress = "results/genomescope2/{sample_id}/progress.txt",
        summary = "results/genomescope2/{sample_id}/summary.txt",
        transformed_linear_plot = "results/genomescope2/{sample_id}/transformed_linear_plot.png",
        transformed_log_plot = "results/genomescope2/{sample_id}/transformed_log_plot.png"
    log:
        out = "logs/genomescope2_{sample_id}.out",
        err = "logs/genomescope2_{sample_id}.err"
    conda:
        "../envs/genomescope2.yml"
    params:
        ploidy = config["ploidy"]
    shell:
        "genomescope2 \
            --input {input} \
            --output $(dirname {output.linear_plot}) \
            --kmer_length 21 \
            --ploidy {params.ploidy} > {log.out} 2> {log.err}"

rule hifiasm:
    input:
        "results/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        bp_hap1_p_ctg_gfa = "results/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.gfa",
        bp_hap1_p_ctg_lowQ_bed = "results/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.lowQ.bed",
        bp_hap1_p_ctg_noseq_gfa = "results/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.noseq.gfa",
        bp_hap2_p_ctg_gfa = "results/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.gfa",
        bp_hap2_p_ctg_lowQ_bed = "results/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.lowQ.bed",
        bp_hap2_p_ctg_noseq_gfa = "results/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.noseq.gfa",
        bp_p_ctg_gfa = "results/hifiasm/{sample_id}.asm.bp.p_ctg.gfa",
        bp_p_ctg_lowQ_bed = "results/hifiasm/{sample_id}.asm.bp.p_ctg.lowQ.bed",
        bp_p_ctg_noseq_gfa = "results/hifiasm/{sample_id}.asm.bp.p_ctg.noseq.gfa",
        bp_p_utg_gfa = "results/hifiasm/{sample_id}.asm.bp.p_utg.gfa",
        bp_p_utg_lowQ_bed = "results/hifiasm/{sample_id}.asm.bp.p_utg.lowQ.bed",
        bp_p_utg_noseq_gfa = "results/hifiasm/{sample_id}.asm.bp.p_utg.noseq.gfa",
        bp_r_utg_gfa = "results/hifiasm/{sample_id}.asm.bp.r_utg.gfa",
        bp_r_utg_lowQ_bed = "results/hifiasm/{sample_id}.asm.bp.r_utg.lowQ.bed",
        bp_r_utg_noseq_gfa = "results/hifiasm/{sample_id}.asm.bp.r_utg.noseq.gfa",
        ec_bin = "results/hifiasm/{sample_id}.asm.ec.bin",
        ovlp_reverse_bin = "results/hifiasm/{sample_id}.asm.ovlp.reverse.bin",
        ovlp_source_bin = "results/hifiasm/{sample_id}.asm.ovlp.source.bin"
    log:
        out = "logs/hifiasm_{sample_id}.out",
        err = "logs/hifiasm_{sample_id}.err"
    conda:
        "../envs/hifiasm.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "hifiasm \
            {input} \
            -o $(dirname {output.bp_p_ctg_gfa})/$(basename {output.bp_p_ctg_gfa} .bp.p_ctg.gfa) \
            -t {threads} > {log.out} 2> {log.err}"

rule convert_gfa_to_fa:
    input:
        "results/hifiasm/{sample_id}.asm.bp.p_ctg.gfa"
    output:
        "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa"
    log:
        "logs/convert_gfa_to_fa_{sample_id}.err"
    conda:
        "../envs/hifiasm.yml"
    shell:
        "awk -f workflow/scripts/convert_gfa_to_fa.awk {input} > {output} 2> {log}"

rule seqkit_stats:
    input:
        "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa"
    output:
        "results/hifiasm/{sample_id}_seqkit_stats.txt"
    log:
        "logs/seqkit_stats_{sample_id}.err"
    conda:
        "../envs/seqkit.yml"
    shell:
        "seqkit stats --all {input} > {output} 2> {log}"

rule download_busco_database:
    output:
        directory("results/busco_downloads/")
    log:
        out = "logs/download_busco_database.out",
        err = "logs/download_busco_database.err"
    conda:
        "../envs/busco.yml"
    params:
        lineage_dataset = config["busco_lineage_dataset"]
    shell:
        """
        (
            busco --download {params.lineage_dataset}
            mv busco_downloads/* {output}
            rm -rf busco_downloads
        ) > {log.out} 2> {log.err}
        """

rule busco_genome_mode:
    input:
        genome = "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa",
        database = "results/busco_downloads"
    output:
        directory("results/busco_genome/{sample_id}")
    log:
        out = "logs/busco_genome_mode_{sample_id}.out",
        err = "logs/busco_genome_mode_{sample_id}.err"
    conda:
        "../envs/busco.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    params:
        lineage_dataset = config["busco_lineage_dataset"]
    shell:
        "busco \
            --in {input.genome} \
            --out_path {output} \
            --mode genome \
            --cpu {threads} \
            --lineage_dataset {params.lineage_dataset} \
            --download_path {input.database} \
            --offline > {log.out} 2> {log.err}"

rule meryl:
    input:
        "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa"
    output:
        "results/merqury/{sample_id}.meryl"
    log:
        out = "logs/meryl_{sample_id}.out",
        err = "logs/meryl_{sample_id}.err"
    conda:
        "../envs/merqury.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "meryl count \
            output {output} \
            k=21 \
            threads {threads} \
            {input} > {log.out} 2> {log.err}"

rule merqury:
    input:
        db = "results/merqury/{sample_id}.meryl",
        assembly = "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa"
    output:
        "results/merqury/{sample_id}.qv"
    log:
        out = "logs/merqury_{sample_id}.out",
        err = "logs/merqury_{sample_id}.err"
    conda:
        "../envs/merqury.yml"
    shell:
        "merqury.sh \
            {input.db} \
            {input.assembly} \
            $(dirname {output})/$(basename {output} .qv) > {log.out} 2> {log.err}"

rule inspector:
    input:
        assembly = "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa",
        reads = "results/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        directory("results/inspector/{sample_id}")
    log:
        out = "logs/inspector_{sample_id}.out",
        err = "logs/inspector_{sample_id}.err"
    conda:
        "../envs/inspector.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "inspector.py \
            --contig {input.assembly} \
            --read {input.reads} \
            --datatype hifi \
            --outpath {output} \
            --thread {threads} > {log.out} 2> {log.err}"

# rule inspector_correct:
#     input:
#         "results/inspector/{sample_id}"
#     output:
#         directory("results/inspector_correct/{sample_id}")
#     log:
#         out = "logs/inspector_correct_{sample_id}.out",
#         err = "logs/inspector_correct_{sample_id}.err"
#     conda:
#         "../envs/inspector.yml"
#     threads:
#         max(1, int(workflow.cores * 0.9))
#     shell:
#         "inspector.py \
#             --inspector {input.assembly} \
#             --datatype pacbio-hifi \
#             --outpath {output} \
#             --thread {threads} > {log.out} 2> {log.err}"
