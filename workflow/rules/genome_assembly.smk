rule bam2fastq:
    input:
        "raw_data/{sample_id}.hifi_reads.bam"
    output:
        "results/hifi_reads/raw_reads/{sample_id}_hifi_reads.fastq.gz"
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
        "results/hifi_reads/raw_reads/{sample_id}_hifi_reads.fastq.gz"
    output:
        reads = "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz",
        report_html = "results/hifi_reads/fastplong/{sample_id}_report.html",
        report_json = "results/hifi_reads/fastplong/{sample_id}_report.json"
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

rule fastk:
    input:
        "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        hist = "results/hifi_reads/smudgeplot/{sample_id}_fastk.hist",
        ktab = "results/hifi_reads/smudgeplot/{sample_id}_fastk.ktab"
    log:
        out = "logs/fastk_{sample_id}.out",
        err = "logs/fastk_{sample_id}.err"
    conda:
        "../envs/smudgeplot.yml"
    threads:
        1
    shell:
        "FastK \
            -v \
            -t4 \
            -k31 \
            -M16 \
            -T{threads} \
            {input} \
            -N$(dirname {output.hist})/$(basename {output.hist} .hist) > {log.out} 2> {log.err}"

rule smudgeplot_hetmers:
    input:
        hist = "results/hifi_reads/smudgeplot/{sample_id}_fastk.hist",
        ktab = "results/hifi_reads/smudgeplot/{sample_id}_fastk.ktab"
    output:
        "results/hifi_reads/smudgeplot/{sample_id}_kmerpairs_text.smu"
    log:
        out = "logs/smudgeplot_hetmers_{sample_id}.out",
        err = "logs/smudgeplot_hetmers_{sample_id}.err"
    conda:
        "../envs/smudgeplot.yml"
    threads:
        1
    shell:
        "smudgeplot.py hetmers \
            -L 12 \
            -t {threads} \
            -o $(dirname {output})/$(basename {output} _text.smu) \
            --verbose \
            $(dirname {input.hist})/$(basename {input.hist} .hist) > {log.out} 2> {log.err}"

rule smudgeplot_all:
    input:
        "results/hifi_reads/smudgeplot/{sample_id}_kmerpairs_text.smu"
    output:
        "results/hifi_reads/smudgeplot/{sample_id}_masked_errors_smu.txt"
    log:
        out = "logs/smudgeplot_all_{sample_id}.out",
        err = "logs/smudgeplot_all_{sample_id}.err"
    conda:
        "../envs/smudgeplot.yml"
    threads:
        1
    shell:
        "smudgeplot.py all \
            -o $(dirname {output})/{wildcards.sample_id} \
            {input} > {log.out} 2> {log.err}"

rule jellyfish_count:
    input:
        "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        "results/hifi_reads/genomescope2/{sample_id}_jellyfish_mer_counts.jf"
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
        "results/hifi_reads/genomescope2/{sample_id}_jellyfish_mer_counts.jf"
    output:
        "results/hifi_reads/genomescope2/{sample_id}_jellyfish.histo"
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
        "results/hifi_reads/genomescope2/{sample_id}_jellyfish.histo"
    output:
        linear_plot = "results/hifi_reads/genomescope2/{sample_id}_linear_plot.png",
        log_plot = "results/hifi_reads/genomescope2/{sample_id}_log_plot.png",
        model = "results/hifi_reads/genomescope2/{sample_id}_model.txt",
        progress = "results/hifi_reads/genomescope2/{sample_id}_progress.txt",
        summary = "results/hifi_reads/genomescope2/{sample_id}_summary.txt",
        transformed_linear_plot = "results/hifi_reads/genomescope2/{sample_id}_transformed_linear_plot.png",
        transformed_log_plot = "results/hifi_reads/genomescope2/{sample_id}_transformed_log_plot.png"
    log:
        out = "logs/genomescope2_{sample_id}.out",
        err = "logs/genomescope2_{sample_id}.err"
    conda:
        "../envs/genomescope2.yml"
    params:
        ploidy = config["ploidy"]
    shell:
        """
        (
            genomescope2 \
                --input {input} \
                --output $(dirname {output.linear_plot}) \
                --kmer_length 21 \
                --ploidy {params.ploidy}
            mv $(dirname {output.linear_plot})/linear_plot.png {output.linear_plot}
            mv $(dirname {output.log_plot})/log_plot.png {output.log_plot}
            mv $(dirname {output.model})/model.txt {output.model}
            mv $(dirname {output.progress})/progress.txt {output.progress}
            mv $(dirname {output.summary})/summary.txt {output.summary}
            mv $(dirname {output.transformed_linear_plot})/transformed_linear_plot.png {output.transformed_linear_plot}
            mv $(dirname {output.transformed_log_plot})/transformed_log_plot.png {output.transformed_log_plot}
        ) > {log.out} 2> {log.err}
        """

rule hifiasm:
    input:
        "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        bp_hap1_p_ctg_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.gfa",
        bp_hap1_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.lowQ.bed",
        bp_hap1_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap1.p_ctg.noseq.gfa",
        bp_hap2_p_ctg_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.gfa",
        bp_hap2_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.lowQ.bed",
        bp_hap2_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.hap2.p_ctg.noseq.gfa",
        bp_p_ctg_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_ctg.gfa",
        bp_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_ctg.lowQ.bed",
        bp_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_ctg.noseq.gfa",
        bp_p_utg_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_utg.gfa",
        bp_p_utg_lowQ_bed = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_utg.lowQ.bed",
        bp_p_utg_noseq_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_utg.noseq.gfa",
        bp_r_utg_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.r_utg.gfa",
        bp_r_utg_lowQ_bed = "results/hifiasm/hifiasm/{sample_id}.asm.bp.r_utg.lowQ.bed",
        bp_r_utg_noseq_gfa = "results/hifiasm/hifiasm/{sample_id}.asm.bp.r_utg.noseq.gfa",
        ec_bin = "results/hifiasm/hifiasm/{sample_id}.asm.ec.bin",
        ovlp_reverse_bin = "results/hifiasm/hifiasm/{sample_id}.asm.ovlp.reverse.bin",
        ovlp_source_bin = "results/hifiasm/hifiasm/{sample_id}.asm.ovlp.source.bin"
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
        "results/hifiasm/hifiasm/{sample_id}.asm.bp.p_ctg.gfa"
    output:
        "results/hifiasm/assembly/{sample_id}.asm.bp.p_ctg.fa"
    log:
        "logs/convert_gfa_to_fa_{sample_id}.err"
    conda:
        "../envs/hifiasm.yml"
    shell:
        "awk -f workflow/scripts/convert_gfa_to_fa.awk {input} > {output} 2> {log}"

rule fcs_gx_get_code:
    output:
        "results/downloads/fcs.py"
    log:
        out = "logs/fcs_gx_get_code.out",
        err = "logs/fcs_gx_get_code.err"
    conda:
        "../envs/fcs_gx.yml" # Dummy
    params:
        url = "https://raw.githubusercontent.com/ncbi/fcs/refs/tags/v0.5.5/dist/fcs.py"
    shell:
        "wget -O {output} {params.url} > {log.out} 2> {log.err}"

rule fcs_gx_get_db:
    input:
        code = "results/downloads/fcs.py"
    output:
        readme = "results/downloads/gxdb/all.README.txt",
        assemblies = "results/downloads/gxdb/all.assemblies.tsv",
        blast_div = "results/downloads/gxdb/all.blast_div.tsv.gz",
        gxi = "results/downloads/gxdb/all.gxi",
        gxs = "results/downloads/gxdb/all.gxs",
        manifest = "results/downloads/gxdb/all.manifest",
        meta = "results/downloads/gxdb/all.meta.jsonl",
        seq_info = "results/downloads/gxdb/all.seq_info.tsv.gz",
        taxa = "results/downloads/gxdb/all.taxa.tsv",
        sif = "results/downloads/fcs_gx_0.5.5.sif"
    log:
        out = "logs/fcs_gx_get_db.out",
        err = "logs/fcs_gx_get_db.err"
    conda:
        "../envs/fcs_gx.yml" # Dummy
    params:
        docker = "docker://ncbi/fcs-gx:0.5.5",
        mft = "https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
    shell:
        """
        (
            singularity pull {output.sif} {params.docker}
            export FCS_DEFAULT_IMAGE={output.sif}
            python3 {input.code} db get \
                --mft {params.mft} \
                --dir $(dirname {output.readme})
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_check_db:
    input:
        code = "results/downloads/fcs.py",
        readme = "results/downloads/gxdb/all.README.txt",
        assemblies = "results/downloads/gxdb/all.assemblies.tsv",
        blast_div = "results/downloads/gxdb/all.blast_div.tsv.gz",
        gxi = "results/downloads/gxdb/all.gxi",
        gxs = "results/downloads/gxdb/all.gxs",
        manifest = "results/downloads/gxdb/all.manifest",
        meta = "results/downloads/gxdb/all.meta.jsonl",
        seq_info = "results/downloads/gxdb/all.seq_info.tsv.gz",
        taxa = "results/downloads/gxdb/all.taxa.tsv",
        sif = "results/downloads/fcs_gx_0.5.5.sif"
    output:
        directory("results/downloads/.gxdb_checked")
    log:
        out = "logs/fcs_gx_check_db.out",
        err = "logs/fcs_gx_check_db.err"
    conda:
        "../envs/fcs_gx.yml" # Dummy
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.sif}
            python3 {input.code} db check \
                --mft {input.manifest} \
                --dir {output}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_screen:
    input:
        code = "results/downloads/fcs.py",
        manifest = "results/downloads/gxdb/all.manifest",
        check = "results/downloads/.gxdb_checked",
        assembly = "results/hifiasm/assembly/{sample_id}.asm.bp.p_ctg.fa",
        sif = "results/downloads/fcs_gx_0.5.5.sif"
    output:
        report = "results/fcs_gx/fcs_gx_screen/{sample_id}.asm.bp.p_ctg." + config["fcs_gx_taxid"] + ".fcs_gx_report.txt",
        taxonomy = "results/fcs_gx/fcs_gx_screen/{sample_id}.asm.bp.p_ctg." + config["fcs_gx_taxid"] + ".taxonomy.rpt"
    log:
        out = "logs/fcs_gx_screen_{sample_id}.out",
        err = "logs/fcs_gx_screen_{sample_id}.err"
    conda:
        "../envs/fcs_gx.yml" # Dummy
    params:
        taxid = config["fcs_gx_taxid"]
    threads:
        48
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.sif}
            python3 {input.code} screen genome \
                --fasta {input.assembly} \
                --out-dir $(dirname {output.report}) \
                --gx-db $(dirname {input.manifest}) \
                --tax-id {params.taxid}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_clean:
    input:
        code = "results/downloads/fcs.py",
        assembly = "results/hifiasm/assembly/{sample_id}.asm.bp.p_ctg.fa",
        screen_report = "results/fcs_gx/fcs_gx_screen/{sample_id}.asm.bp.p_ctg." + config["fcs_gx_taxid"] + ".fcs_gx_report.txt",
        sif = "results/downloads/fcs_gx_0.5.5.sif"
    output:
        clean = "results/fcs_gx/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa",
        contam = "results/fcs_gx/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.contam.fa"
    log:
        out = "logs/fcs_gx_clean_{sample_id}.out",
        err = "logs/fcs_gx_clean_{sample_id}.err"
    conda:
        "../envs/fcs_gx.yml" # Dummy
    params:
        taxid = config["fcs_gx_taxid"]
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.sif}
            cat {input.assembly} | \
            python3 {input.code} clean genome \
                --action-report {input.screen_report} \
                --output {output.clean} \
                --contam-fasta-out {output.contam}
        ) > {log.out} 2> {log.err}
        """

rule copy_fcs_gx_clean:
    input:
        "results/fcs_gx/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
    output:
        "results/fcs_gx/assembly/{sample_id}.asm.bp.p_ctg.fa"
    log:
        out = "logs/copy_fcs_gx_clean_{sample_id}.out",
        err = "logs/copy_fcs_gx_clean_{sample_id}.err"
    conda:
        "../envs/fcs_gx.yml"
    shell:
        "cp {input} {output} > {log.out} 2> {log.err}"

rule seqkit_stats:
    input:
        "results/{assembly}/assembly/{sample_id}.asm.bp.p_ctg.fa"
    output:
        "results/{assembly}/seqkit/{sample_id}_seqkit_stats.txt"
    log:
        "logs/seqkit_stats_{assembly}_{sample_id}.err"
    conda:
        "../envs/seqkit.yml"
    shell:
        "seqkit stats --all {input} > {output} 2> {log}"

rule download_busco_database:
    output:
        directory("results/downloads/busco_downloads")
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
            mkdir -p {output}
            mv busco_downloads/* {output}
            rm -rf busco_downloads
        ) > {log.out} 2> {log.err}
        """

rule busco_genome_mode:
    input:
        genome = "results/{assembly}/assembly/{sample_id}.asm.bp.p_ctg.fa",
        database = "results/downloads/busco_downloads"
    output:
        directory("results/{assembly}/busco_genome/BUSCO_{sample_id}.asm.bp.p_ctg.fa")
    log:
        out = "logs/busco_genome_mode_{assembly}_{sample_id}.out",
        err = "logs/busco_genome_mode_{assembly}_{sample_id}.err"
    conda:
        "../envs/busco.yml"
    threads:
        max(1, int(workflow.cores * 0.8))
    params:
        lineage_dataset = config["busco_lineage_dataset"]
    shell:
        "busco \
            --in {input.genome} \
            --out_path $(dirname {output}) \
            --mode genome \
            --cpu {threads} \
            --lineage_dataset {params.lineage_dataset} \
            --download_path {input.database} \
            --offline > {log.out} 2> {log.err}"

rule meryl:
    input:
        "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        directory("results/hifi_reads/meryl/{sample_id}")
    log:
        out = "logs/meryl_{sample_id}.out",
        err = "logs/meryl_{sample_id}.err"
    conda:
        "../envs/merqury.yml"
    threads:
        1
    shell:
        "meryl count \
            output {output} \
            k=21 \
            threads={threads} \
            {input} > {log.out} 2> {log.err}"

rule merqury:
    input:
        db = "results/hifi_reads/meryl/{sample_id}",
        assembly = "results/{assembly}/assembly/{sample_id}.asm.bp.p_ctg.fa"
    output:
        "results/{assembly}/merqury/{sample_id}.merqury.qv"
    log:
        out = "logs/merqury_{assembly}_{sample_id}.out",
        err = "logs/merqury_{assembly}_{sample_id}.err"
    conda:
        "../envs/merqury.yml"
    shell:
        """
        (
            cd $(dirname {output})
            merqury.sh \
                ../../../{input.db} \
                ../../../{input.assembly} \
                $(basename {output} .qv)
            cd ../../..
        ) > {log.out} 2> {log.err}
        """

rule inspector:
    input:
        assembly = "results/{assembly}/assembly/{sample_id}.asm.bp.p_ctg.fa",
        reads = "results/hifi_reads/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        directory("results/{assembly}/inspector/{sample_id}")
    log:
        out = "logs/inspector_{assembly}_{sample_id}.out",
        err = "logs/inspector_{assembly}_{sample_id}.err"
    conda:
        "../envs/inspector.yml"
    threads:
        max(1, int(workflow.cores * 0.8))
    shell:
        "inspector.py \
            --contig {input.assembly} \
            --read {input.reads} \
            --datatype hifi \
            --outpath {output} \
            --thread {threads} > {log.out} 2> {log.err}"

# rule inspector_correct:
#     input:
#         "results/{assembly}/inspector/{sample_id}"
#     output:
#         directory("results/{assembly}/inspector_correct/{sample_id}")
#     log:
#         out = "logs/inspector_correct_{assembly}_{sample_id}.out",
#         err = "logs/inspector_correct_{assembly}_{sample_id}.err"
#     conda:
#         "../envs/inspector.yml"
#     threads:
#         max(1, int(workflow.cores * 0.9))
#     shell:
#         "inspector.py \
#             --inspector {input} \
#             --datatype pacbio-hifi \
#             --outpath {output} \
#             --thread {threads} > {log.out} 2> {log.err}"

rule gt_suffixerator:
    input:
        "results/{assembly}/assembly/{sample_id}.asm.bp.p_ctg.fa"
    output:
        assembly = "results/{assembly}/lai/gt_suffixerator/{sample_id}.fa",
        des = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.des",
        esq = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.esq",
        lcp = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.lcp",
        llv = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.llv",
        md5 = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.md5",
        prj = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.prj",
        sds = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.sds",
        ssp = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.ssp",
        suf = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.suf"
    log:
        out = "logs/gt_suffixerator_{assembly}_{sample_id}.out",
        err = "logs/gt_suffixerator_{assembly}_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        """
        (
            cp {input} {output.assembly}
            gt suffixerator \
                -db {output.assembly} \
                -indexname $(dirname {output.des})/$(basename {output.des} .des) \
                -tis \
                -suf \
                -lcp \
                -des \
                -ssp \
                -sds \
                -dna
        ) > {log.out} 2> {log.err}
        """

rule gt_ltrharvest:
    input:
        assembly = "results/{assembly}/lai/gt_suffixerator/{sample_id}.fa",
        des = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.des",
        esq = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.esq",
        lcp = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.lcp",
        llv = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.llv",
        md5 = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.md5",
        prj = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.prj",
        sds = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.sds",
        ssp = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.ssp",
        suf = "results/{assembly}/lai/gt_suffixerator/{sample_id}_index.suf"
    output:
        "results/{assembly}/lai/gt_ltrharvest/{sample_id}.fa.harvest.scn"
    log:
        err = "logs/gt_ltrharvest_{assembly}_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        "gt ltrharvest \
            -index $(dirname {input.des})/$(basename {input.des} .des) \
            -minlenltr 100 \
            -maxlenltr 7000 \
            -mintsd 4 \
            -maxtsd 6 \
            -motif TGCA \
            -motifmis 1 \
            -similar 85 \
            -vic 10 \
            -seed 20 \
            -seqids yes > {output} 2> {log.err}"

rule ltr_finder_parallel:
    input:
        "results/{assembly}/lai/gt_suffixerator/{sample_id}.fa"
    output:
        "results/{assembly}/lai/ltr_finder_parallel/{sample_id}.fa.finder.combine.scn"
    log:
        out = "logs/ltr_finder_parallel_{assembly}_{sample_id}.out",
        err = "logs/ltr_finder_parallel_{assembly}_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    threads:
        max(1, int(workflow.cores * 0.1) - 1)
    shell:
        """
        (
            cd $(dirname {output})
            LTR_FINDER_parallel \
                -seq ../../../../{input} \
                -threads {threads} \
                -harvest_out \
                -size 1000000 \
                -time 300
            cd ../../../../
        ) > {log.out} 2> {log.err}
        """

rule merge_ltrharvest_and_ltrfinder:
    input:
        ltrharvest = "results/{assembly}/lai/gt_ltrharvest/{sample_id}.fa.harvest.scn",
        ltrfinder = "results/{assembly}/lai/ltr_finder_parallel/{sample_id}.fa.finder.combine.scn"
    output:
        "results/{assembly}/lai/lib_ltrharvest_ltrfinder/{sample_id}.fa.rawLTR.scn"
    log:
        err = "logs/merge_ltrharvest_and_ltrfinder_{assembly}_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        "cat {input.ltrharvest} {input.ltrfinder} > {output} 2> {log.err}"

rule ltr_retriever:
    input:
        assembly = "results/{assembly}/lai/gt_suffixerator/{sample_id}.fa",
        inharvest = "results/{assembly}/lai/lib_ltrharvest_ltrfinder/{sample_id}.fa.rawLTR.scn"
    output:
        "results/{assembly}/lai/ltr_retriever/{sample_id}.fa.out.LAI"
    log:
        out = "logs/ltr_retriever_{assembly}_{sample_id}.out",
        err = "logs/ltr_retriever_{assembly}_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    threads:
        max(1, int(workflow.cores * 0.1) - 1)
    shell:
        """
        (
            cd $(dirname {output})
            LTR_retriever \
                -genome ../../../../{input.assembly} \
                -inharvest ../../../../{input.inharvest} \
                -threads {threads}
            cd ../../../../
        ) > {log.out} 2> {log.err}
        """
