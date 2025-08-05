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

rule fastk:
    input:
        "results/fastplong/{sample_id}_hifi_reads_curated.fastq.gz"
    output:
        hist = "results/fastk/{sample_id}.hist",
        ktab = "results/fastk/{sample_id}.ktab"
    log:
        out = "logs/fastk_{sample_id}.out",
        err = "logs/fastk_{sample_id}.err"
    conda:
        "../envs/fastk.yml"
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
        hist = "results/fastk/{sample_id}.hist",
        ktab = "results/fastk/{sample_id}.ktab"
    output:
        "results/fastk/{sample_id}_kmerpairs_text.smu"
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
            {input} > {log.out} 2> {log.err}"

rule smudgeplot_all:
    input:
        "results/fastk/{sample_id}_kmerpairs_text.smu"
    output:
        directory("results/smudgeplot/{sample_id}")
    log:
        out = "logs/smudgeplot_all_{sample_id}.out",
        err = "logs/smudgeplot_all_{sample_id}.err"
    conda:
        "../envs/smudgeplot.yml"
    threads:
        1
    shell:
        "smudgeplot.py all \
            -o {output}/smudgeplot \
            {input} > {log.out} 2> {log.err}"

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

rule fcs_gx_get_code:
    output:
        "results/fcs.py"
    log:
        out = "logs/fcs_gx_get_code.out",
        err = "logs/fcs_gx_get_code.err"
    conda:
        "../envs/fcs-gx.yml" # Dummy
    shell:
        "wget -O {output} https://raw.githubusercontent.com/ncbi/fcs/refs/tags/v0.5.5/dist/fcs.py > {log.out} 2> {log.err}"

rule fcs_gx_get_db:
    input:
        code = "results/fcs.py"
    output:
        readme = "results/fcs_gx_db/gxdb/all.README.txt",
        assemblies = "results/fcs_gx_db/gxdb/all.assemblies.tsv",
        blast_div = "results/fcs_gx_db/gxdb/all.blast_div.tsv.gz",
        gxi = "results/fcs_gx_db/gxdb/all.gxi",
        gxs = "results/fcs_gx_db/gxdb/all.gxs",
        manifest = "results/fcs_gx_db/gxdb/all.manifest",
        meta = "results/fcs_gx_db/gxdb/all.meta.jsonl",
        seq_info = "results/fcs_gx_db/gxdb/all.seq_info.tsv.gz",
        taxa = "results/fcs_gx_db/gxdb/all.taxa.tsv",
        sif = "results/fcs-gx_0.5.5.sif"
    log:
        out = "logs/fcs_gx_get_db.out",
        err = "logs/fcs_gx_get_db.err"
    conda:
        "../envs/fcs-gx.yml" # Dummy
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
                --dir $(dirname {output.readme}) > {log.out} 2> {log.err}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_check_db:
    input:
        code = "results/fcs.py",
        readme = "results/fcs_gx_db/gxdb/all.README.txt",
        assemblies = "results/fcs_gx_db/gxdb/all.assemblies.tsv",
        blast_div = "results/fcs_gx_db/gxdb/all.blast_div.tsv.gz",
        gxi = "results/fcs_gx_db/gxdb/all.gxi",
        gxs = "results/fcs_gx_db/gxdb/all.gxs",
        manifest = "results/fcs_gx_db/gxdb/all.manifest",
        meta = "results/fcs_gx_db/gxdb/all.meta.jsonl",
        seq_info = "results/fcs_gx_db/gxdb/all.seq_info.tsv.gz",
        taxa = "results/fcs_gx_db/gxdb/all.taxa.tsv",
        sif = "results/fcs-gx_0.5.5.sif"
    output:
        directory("results/fcs_gx_db/check")
    log:
        out = "logs/fcs_gx_check_db.out",
        err = "logs/fcs_gx_check_db.err"
    conda:
        "../envs/fcs-gx.yml" # Dummy
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
        code = "results/fcs.py",
        manifest = "results/fcs_gx_db/gxdb/all.manifest",
        check = "results/fcs_gx_db/check",
        assembly = "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa",
        sif = "results/fcs-gx_0.5.5.sif"
    output:
        directory("results/fcs_gx_screen/{sample_id}")
    log:
        out = "logs/fcs_gx_screen_{sample_id}.out",
        err = "logs/fcs_gx_screen_{sample_id}.err"
    conda:
        "../envs/fcs-gx.yml" # Dummy
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
                --out-dir {output} \
                --gx-db $(dirname {input.manifest}) \
                --tax-id {params.taxid}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_clean:
    input:
        code = "results/fcs.py",
        assembly = "results/hifiasm/{sample_id}.asm.bp.p_ctg.fa",
        screen = "results/fcs_gx_screen/{sample_id}",
        sif = "results/fcs-gx_0.5.5.sif"
    output:
        clean = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa",
        contam = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.contam.fa"
    log:
        out = "logs/fcs_gx_clean_{sample_id}.out",
        err = "logs/fcs_gx_clean_{sample_id}.err"
    conda:
        "../envs/fcs-gx.yml" # Dummy
    params:
        taxid = config["fcs_gx_taxid"]
    threads:
        48
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.sif}
            zcat {input.assembly} | \
            python3 {input.code} clean genome \
                --action-report {input.screen}/{wildcards.sample_id}.asm.bp.p_ctg.fa.{params.taxid}.fcs_gx_report.txt \
                --output {output.clean} \
                --contam-fasta-out {output.contam}
        ) > {log.out} 2> {log.err}
        """

rule seqkit_stats:
    input:
        "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
    output:
        "results/fcs_gx_clean/{sample_id}_seqkit_stats.txt"
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
            mkdir -p {output}
            mv busco_downloads/* {output}
            rm -rf busco_downloads
        ) > {log.out} 2> {log.err}
        """

rule busco_genome_mode:
    input:
        genome = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa",
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
        "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
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
        assembly = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
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
        assembly = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa",
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

rule gt_suffixerator:
    input:
        "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
    output:
        assembly = "results/lai/{sample_id}.fa",
        index = "results/lai/{sample_id}_index"
    log:
        out = "logs/gt_suffixerator_{sample_id}.out",
        err = "logs/gt_suffixerator_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        """
        (
            cp {input} {output.assembly}
            gt suffixerator \
                -db {output.assembly} \
                -indexname {output.index} \
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
        "results/lai/{sample_id}_index"
    output:
        "results/lai/{sample_id}.fa.harvest.scn"
    log:
        err = "logs/gt_ltrharvest_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        "gt ltrharvest \
            -index {input} \
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
        "results/lai/{sample_id}.fa"
    output:
        "results/lai/{sample_id}.fa.finder.combine.scn"
    log:
        out = "logs/ltr_finder_parallel_{sample_id}.out",
        err = "logs/ltr_finder_parallel_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "LTR_FINDER_parallel \
            -seq {input} \
            -threads {threads} \
            -harvest_out \
            -size 1000000 \
            -time 300 > {log.out} 2> {log.err}"

rule merge_ltrharvest_and_ltrfinder:
    input:
        ltrharvest = "results/lai/{sample_id}.fa.harvest.scn",
        ltrfinder = "results/lai/{sample_id}.fa.finder.combine.scn"
    output:
        "results/lai/{sample_id}.fa.rawLTR.scn"
    log:
        err = "logs/merge_ltrharvest_and_ltrfinder_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    shell:
        "cat {input.ltrharvest} {input.ltrfinder} > {output} 2> {log.err}"

rule ltr_retriever:
    input:
        assembly = "results/lai/{sample_id}.fa",
        inharvest = "results/lai/{sample_id}.fa.rawLTR.scn"
    output:
        "results/lai/{sample_id}.fa.pass.list"
    log:
        out = "logs/ltr_retriever_{sample_id}.out",
        err = "logs/ltr_retriever_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "LTR_retriever \
            -genome {input.assembly} \
            -inharvest {input.inharvest} \
            -threads {threads} > {log.out} 2> {log.err}"

rule lai:
    input:
        assembly = "results/lai/{sample_id}.fa",
        intact = "results/lai/{sample_id}.fa.pass.list"
    output:
        "results/lai/{sample_id}.fa.out.LAI"
    log:
        out = "logs/lai_{sample_id}.out",
        err = "logs/lai_{sample_id}.err"
    conda:
        "../envs/lai.yml"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "LAI \
            -genome {input.assembly} \
            -intact {input.intact} \
            -all {input.assembly}.out \
            -t {threads} > {log.out} 2> {log.err}"
