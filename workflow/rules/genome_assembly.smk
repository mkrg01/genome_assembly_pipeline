import os
import glob
hifi_samples = glob.glob("raw_data/*.hifi_reads.bam")
if len(hifi_samples) == 0:
    raise RuntimeError("No HiFi read files matching 'raw_data/*.hifi_reads.bam' were found. Please check the input directory.")
hifi_sample_ids = sorted({os.path.basename(f).replace(".hifi_reads.bam", "") for f in hifi_samples})

wildcard_constraints:
    assembly_name = config["assembly_name"]

rule bam2fastq:
    input:
        read = "raw_data/{hifi_sample_id}.hifi_reads.bam",
        index = "raw_data/{hifi_sample_id}.hifi_reads.bam.pbi"
    output:
        "results/hifi_reads/raw_reads/{hifi_sample_id}_hifi_reads.fastq.gz"
    log:
        out = "logs/bam2fastq_{hifi_sample_id}.out",
        err = "logs/bam2fastq_{hifi_sample_id}.err"
    conda:
        "../envs/bam2fastq.yml"
    threads:
        max(1, int(workflow.cores * 0.1))
    shell:
        "bam2fastq \
            {input.read} \
            --output $(dirname {output})/$(basename {output} .fastq.gz) \
            --num-threads {threads} > {log.out} 2> {log.err}"

rule fastplong:
    input:
        "results/hifi_reads/raw_reads/{hifi_sample_id}_hifi_reads.fastq.gz"
    output:
        reads = "results/hifi_reads/fastplong/{hifi_sample_id}_hifi_reads_curated.fastq.gz",
        report_html = "results/hifi_reads/fastplong/{hifi_sample_id}_report.html",
        report_json = "results/hifi_reads/fastplong/{hifi_sample_id}_report.json"
    log:
        out = "logs/fastplong_{hifi_sample_id}.out",
        err = "logs/fastplong_{hifi_sample_id}.err"
    conda:
        "../envs/fastplong.yml"
    threads:
        max(1, int(workflow.cores * 0.1))
    shell:
        "fastplong \
            --in {input} \
            --out {output.reads} \
            --html {output.report_html} \
            --json {output.report_json} \
            --thread {threads} > {log.out} 2> {log.err}"

rule merge_or_copy_hifi_reads:
    input:
        expand("results/hifi_reads/fastplong/{hifi_sample_id}_hifi_reads_curated.fastq.gz", hifi_sample_id=hifi_sample_ids)
    output:
        "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    log:
        out = "logs/fastplong_{assembly_name}.out",
        err = "logs/fastplong_{assembly_name}.err"
    conda:
        "../envs/fastplong.yml"
    shell:
        "cat {input} > {output} 2> {log.err}"

rule fastk:
    input:
        "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        hist = "results/hifi_reads/smudgeplot/{assembly_name}_fastk.hist",
        ktab = "results/hifi_reads/smudgeplot/{assembly_name}_fastk.ktab"
    log:
        out = "logs/fastk_{assembly_name}.out",
        err = "logs/fastk_{assembly_name}.err"
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
        hist = "results/hifi_reads/smudgeplot/{assembly_name}_fastk.hist",
        ktab = "results/hifi_reads/smudgeplot/{assembly_name}_fastk.ktab"
    output:
        "results/hifi_reads/smudgeplot/{assembly_name}_kmerpairs_text.smu"
    log:
        out = "logs/smudgeplot_hetmers_{assembly_name}.out",
        err = "logs/smudgeplot_hetmers_{assembly_name}.err"
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
        "results/hifi_reads/smudgeplot/{assembly_name}_kmerpairs_text.smu"
    output:
        "results/hifi_reads/smudgeplot/{assembly_name}_masked_errors_smu.txt"
    log:
        out = "logs/smudgeplot_all_{assembly_name}.out",
        err = "logs/smudgeplot_all_{assembly_name}.err"
    conda:
        "../envs/smudgeplot.yml"
    threads:
        1
    shell:
        "smudgeplot.py all \
            -o $(dirname {output})/{wildcards.assembly_name} \
            {input} > {log.out} 2> {log.err}"

rule jellyfish_count:
    input:
        "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        "results/hifi_reads/genomescope2/{assembly_name}_jellyfish_mer_counts.jf"
    log:
        out = "logs/jellyfish_count_{assembly_name}.out",
        err = "logs/jellyfish_count_{assembly_name}.err"
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
        "results/hifi_reads/genomescope2/{assembly_name}_jellyfish_mer_counts.jf"
    output:
        "results/hifi_reads/genomescope2/{assembly_name}_jellyfish.histo"
    log:
        out = "logs/jellyfish_histo_{assembly_name}.out",
        err = "logs/jellyfish_histo_{assembly_name}.err"
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
        "results/hifi_reads/genomescope2/{assembly_name}_jellyfish.histo"
    output:
        linear_plot = "results/hifi_reads/genomescope2/{assembly_name}_linear_plot.png",
        log_plot = "results/hifi_reads/genomescope2/{assembly_name}_log_plot.png",
        model = "results/hifi_reads/genomescope2/{assembly_name}_model.txt",
        progress = "results/hifi_reads/genomescope2/{assembly_name}_progress.txt",
        summary = "results/hifi_reads/genomescope2/{assembly_name}_summary.txt",
        transformed_linear_plot = "results/hifi_reads/genomescope2/{assembly_name}_transformed_linear_plot.png",
        transformed_log_plot = "results/hifi_reads/genomescope2/{assembly_name}_transformed_log_plot.png"
    log:
        out = "logs/genomescope2_{assembly_name}.out",
        err = "logs/genomescope2_{assembly_name}.err"
    conda:
        "../envs/genomescope2.yml"
    shell:
        """
        (
            genomescope2 \
                --input {input} \
                --output $(dirname {output.linear_plot}) \
                --kmer_length 21
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
        "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        bp_hap1_p_ctg_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap1.p_ctg.gfa",
        bp_hap1_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap1.p_ctg.lowQ.bed",
        bp_hap1_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap1.p_ctg.noseq.gfa",
        bp_hap2_p_ctg_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap2.p_ctg.gfa",
        bp_hap2_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap2.p_ctg.lowQ.bed",
        bp_hap2_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap2.p_ctg.noseq.gfa",
        bp_p_ctg_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_ctg.gfa",
        bp_p_ctg_lowQ_bed = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_ctg.lowQ.bed",
        bp_p_ctg_noseq_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_ctg.noseq.gfa",
        bp_p_utg_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_utg.gfa",
        bp_p_utg_lowQ_bed = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_utg.lowQ.bed",
        bp_p_utg_noseq_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_utg.noseq.gfa",
        bp_r_utg_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.r_utg.gfa",
        bp_r_utg_lowQ_bed = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.r_utg.lowQ.bed",
        bp_r_utg_noseq_gfa = "results/hifiasm/hifiasm/{assembly_name}.asm.bp.r_utg.noseq.gfa",
        ec_bin = "results/hifiasm/hifiasm/{assembly_name}.asm.ec.bin",
        ovlp_reverse_bin = "results/hifiasm/hifiasm/{assembly_name}.asm.ovlp.reverse.bin",
        ovlp_source_bin = "results/hifiasm/hifiasm/{assembly_name}.asm.ovlp.source.bin"
    log:
        out = "logs/hifiasm_{assembly_name}.out",
        err = "logs/hifiasm_{assembly_name}.err"
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
        "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_ctg.gfa"
    output:
        "results/hifiasm/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    log:
        "logs/convert_gfa_to_fa_{assembly_name}.err"
    conda:
        "../envs/hifiasm.yml"
    shell:
        "awk -f workflow/scripts/convert_gfa_to_fa.awk {input} > {output} 2> {log}"

rule fcs_get_code:
    output:
        fcs_py = "results/downloads/fcs/fcs.py",
        run_fcsadaptor_sh = "results/downloads/fcs/run_fcsadaptor.sh"
    log:
        out = "logs/fcs_get_code.out",
        err = "logs/fcs_get_code.err"
    conda:
        "../envs/fcs.yml"
    params:
        fcs_py_url = "https://raw.githubusercontent.com/ncbi/fcs/refs/tags/v0.5.5/dist/fcs.py",
        run_fcsadaptor_sh_url = "https://raw.githubusercontent.com/ncbi/fcs/refs/tags/v0.5.5/dist/run_fcsadaptor.sh"
    shell:
        """
        (
            wget -O {output.fcs_py} {params.fcs_py_url}
            wget -O {output.run_fcsadaptor_sh} {params.run_fcsadaptor_sh_url}
            chmod 755 {output.run_fcsadaptor_sh}
        ) > {log.out} 2> {log.err}
        """

rule fcs_pull_image:
    output:
        fcs_adaptor_sif = "results/downloads/fcs/fcs_adaptor_0.5.5.sif",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    log:
        out = "logs/fcs_pull_image.out",
        err = "logs/fcs_pull_image.err"
    conda:
        "../envs/fcs.yml"
    params:
        fcs_adaptor_docker = "docker://ncbi/fcs-adaptor:0.5.5",
        fcs_gx_docker = "docker://ncbi/fcs-gx:0.5.5"
    shell:
        """
        (
            singularity pull {output.fcs_adaptor_sif} {params.fcs_adaptor_docker}
            singularity pull {output.fcs_gx_sif} {params.fcs_gx_docker}
        ) > {log.out} 2> {log.err}
        """

rule fcs_adaptor_screen:
    input:
        run_fcsadaptor_sh = "results/downloads/fcs/run_fcsadaptor.sh",
        assembly = "results/hifiasm/assembly/{assembly_name}.asm.bp.p_ctg.fa",
        fcs_adaptor_sif = "results/downloads/fcs/fcs_adaptor_0.5.5.sif"
    output:
        assembly = "results/fcs/fcs_adaptor_screen/cleaned_sequences/{assembly_name}.asm.bp.p_ctg.fa",
        combined_calls_jsonl = "results/fcs/fcs_adaptor_screen/{assembly_name}_combined.calls.jsonl",
        fcs_log = "results/fcs/fcs_adaptor_screen/{assembly_name}_fcs.log",
        fcs_adaptor_log = "results/fcs/fcs_adaptor_screen/{assembly_name}_fcs_adaptor.log",
        fcs_adaptor_report_txt = "results/fcs/fcs_adaptor_screen/{assembly_name}_fcs_adaptor_report.txt",
        logs_jsonl = "results/fcs/fcs_adaptor_screen/{assembly_name}_logs.jsonl",
        pipeline_args_yaml = "results/fcs/fcs_adaptor_screen/{assembly_name}_pipeline_args.yaml",
        skipped_trims_jsonl = "results/fcs/fcs_adaptor_screen/{assembly_name}_skipped_trims.jsonl",
        validate_fasta_txt = "results/fcs/fcs_adaptor_screen/{assembly_name}_validate_fasta.txt"
    log:
        out = "logs/fcs_adaptor_screen_{assembly_name}.out",
        err = "logs/fcs_adaptor_screen_{assembly_name}.err"
    conda:
        "../envs/fcs.yml"
    shell:
        """
        (
            mkdir -p $(dirname {output.fcs_adaptor_report_txt})
            {input.run_fcsadaptor_sh} \
                --fasta-input {input.assembly} \
                --output-dir $(dirname {output.fcs_adaptor_report_txt}) \
                --euk \
                --container-engine singularity \
                --image {input.fcs_adaptor_sif}
            mv $(dirname {output.combined_calls_jsonl})/combined.calls.jsonl {output.combined_calls_jsonl}
            mv $(dirname {output.fcs_log})/fcs.log {output.fcs_log}
            mv $(dirname {output.fcs_adaptor_log})/fcs_adaptor.log {output.fcs_adaptor_log}
            mv $(dirname {output.fcs_adaptor_report_txt})/fcs_adaptor_report.txt {output.fcs_adaptor_report_txt}
            mv $(dirname {output.logs_jsonl})/logs.jsonl {output.logs_jsonl}
            mv $(dirname {output.pipeline_args_yaml})/pipeline_args.yaml {output.pipeline_args_yaml}
            mv $(dirname {output.skipped_trims_jsonl})/skipped_trims.jsonl {output.skipped_trims_jsonl}
            mv $(dirname {output.validate_fasta_txt})/validate_fasta.txt {output.validate_fasta_txt}
        ) > {log.out} 2> {log.err}
        """

rule fcs_adaptor_clean:
    input:
        fcs_py = "results/downloads/fcs/fcs.py",
        assembly = "results/hifiasm/assembly/{assembly_name}.asm.bp.p_ctg.fa",
        fcs_adaptor_report_txt = "results/fcs/fcs_adaptor_screen/{assembly_name}_fcs_adaptor_report.txt",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    output:
        clean = "results/fcs/fcs_adaptor_clean/{assembly_name}.asm.bp.p_ctg.clean.fa",
        contam = "results/fcs/fcs_adaptor_clean/{assembly_name}.asm.bp.p_ctg.contam.fa"
    log:
        out = "logs/fcs_adaptor_clean_{assembly_name}.out",
        err = "logs/fcs_adaptor_clean_{assembly_name}.err"
    conda:
        "../envs/fcs.yml"
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.fcs_gx_sif}
            cat {input.assembly} | \
            python3 {input.fcs_py} clean genome \
                --action-report {input.fcs_adaptor_report_txt} \
                --output {output.clean} \
                --contam-fasta-out {output.contam} \
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_get_db:
    input:
        code = "results/downloads/fcs/fcs.py",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    output:
        readme = "results/downloads/fcs/gxdb/all.README.txt",
        assemblies = "results/downloads/fcs/gxdb/all.assemblies.tsv",
        blast_div = "results/downloads/fcs/gxdb/all.blast_div.tsv.gz",
        gxi = "results/downloads/fcs/gxdb/all.gxi",
        gxs = "results/downloads/fcs/gxdb/all.gxs",
        manifest = "results/downloads/fcs/gxdb/all.manifest",
        meta = "results/downloads/fcs/gxdb/all.meta.jsonl",
        seq_info = "results/downloads/fcs/gxdb/all.seq_info.tsv.gz",
        taxa = "results/downloads/fcs/gxdb/all.taxa.tsv"
    log:
        out = "logs/fcs_gx_get_db.out",
        err = "logs/fcs_gx_get_db.err"
    conda:
        "../envs/fcs.yml"
    params:
        mft = "https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.fcs_gx_sif}
            python3 {input.code} db get \
                --mft {params.mft} \
                --dir $(dirname {output.readme})
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_check_db:
    input:
        code = "results/downloads/fcs/fcs.py",
        readme = "results/downloads/fcs/gxdb/all.README.txt",
        assemblies = "results/downloads/fcs/gxdb/all.assemblies.tsv",
        blast_div = "results/downloads/fcs/gxdb/all.blast_div.tsv.gz",
        gxi = "results/downloads/fcs/gxdb/all.gxi",
        gxs = "results/downloads/fcs/gxdb/all.gxs",
        manifest = "results/downloads/fcs/gxdb/all.manifest",
        meta = "results/downloads/fcs/gxdb/all.meta.jsonl",
        seq_info = "results/downloads/fcs/gxdb/all.seq_info.tsv.gz",
        taxa = "results/downloads/fcs/gxdb/all.taxa.tsv",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    output:
        directory("results/downloads/fcs/.gxdb_checked")
    log:
        out = "logs/fcs_gx_check_db.out",
        err = "logs/fcs_gx_check_db.err"
    conda:
        "../envs/fcs.yml"
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.fcs_gx_sif}
            python3 {input.code} db check \
                --mft {input.manifest} \
                --dir {output}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_screen:
    input:
        code = "results/downloads/fcs/fcs.py",
        manifest = "results/downloads/fcs/gxdb/all.manifest",
        check = "results/downloads/fcs/.gxdb_checked",
        assembly = "results/fcs/fcs_adaptor_clean/{assembly_name}.asm.bp.p_ctg.clean.fa",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    output:
        report = f"results/fcs/fcs_gx_screen/{{assembly_name}}.asm.bp.p_ctg.clean.{config["fcs_gx_taxid"]}.fcs_gx_report.txt",
        taxonomy = f"results/fcs/fcs_gx_screen/{{assembly_name}}.asm.bp.p_ctg.clean.{config["fcs_gx_taxid"]}.taxonomy.rpt"
    log:
        out = "logs/fcs_gx_screen_{assembly_name}.out",
        err = "logs/fcs_gx_screen_{assembly_name}.err"
    conda:
        "../envs/fcs.yml"
    params:
        taxid = config["fcs_gx_taxid"]
    threads:
        48
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.fcs_gx_sif}
            python3 {input.code} screen genome \
                --fasta {input.assembly} \
                --out-dir $(dirname {output.report}) \
                --gx-db $(dirname {input.manifest}) \
                --tax-id {params.taxid}
        ) > {log.out} 2> {log.err}
        """

rule fcs_gx_clean:
    input:
        code = "results/downloads/fcs/fcs.py",
        assembly = "results/fcs/fcs_adaptor_clean/{assembly_name}.asm.bp.p_ctg.clean.fa",
        screen_report = f"results/fcs/fcs_gx_screen/{{assembly_name}}.asm.bp.p_ctg.clean.{config["fcs_gx_taxid"]}.fcs_gx_report.txt",
        fcs_gx_sif = "results/downloads/fcs/fcs_gx_0.5.5.sif"
    output:
        clean = "results/fcs/fcs_gx_clean/{assembly_name}.asm.bp.p_ctg.clean.fa",
        contam = "results/fcs/fcs_gx_clean/{assembly_name}.asm.bp.p_ctg.contam.fa"
    log:
        out = "logs/fcs_gx_clean_{assembly_name}.out",
        err = "logs/fcs_gx_clean_{assembly_name}.err"
    conda:
        "../envs/fcs.yml"
    params:
        taxid = config["fcs_gx_taxid"]
    shell:
        """
        (
            export FCS_DEFAULT_IMAGE={input.fcs_gx_sif}
            cat {input.assembly} | \
            python3 {input.code} clean genome \
                --action-report {input.screen_report} \
                --output {output.clean} \
                --contam-fasta-out {output.contam}
        ) > {log.out} 2> {log.err}
        """

rule copy_fcs_gx_clean:
    input:
        "results/fcs/fcs_gx_clean/{assembly_name}.asm.bp.p_ctg.clean.fa"
    output:
        "results/fcs/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    log:
        out = "logs/copy_fcs_gx_clean_{assembly_name}.out",
        err = "logs/copy_fcs_gx_clean_{assembly_name}.err"
    conda:
        "../envs/fcs.yml"
    shell:
        "cp {input} {output} > {log.out} 2> {log.err}"

rule seqkit_stats:
    input:
        "results/{assembly}/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        "results/{assembly}/seqkit/{assembly_name}_seqkit_stats.txt"
    log:
        "logs/seqkit_stats_{assembly}_{assembly_name}.err"
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
        genome = "results/{assembly}/assembly/{assembly_name}.asm.bp.p_ctg.fa",
        database = "results/downloads/busco_downloads"
    output:
        directory("results/{assembly}/busco_genome/BUSCO_{assembly_name}.asm.bp.p_ctg.fa")
    log:
        out = "logs/busco_genome_mode_{assembly}_{assembly_name}.out",
        err = "logs/busco_genome_mode_{assembly}_{assembly_name}.err"
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
        "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        directory("results/hifi_reads/meryl/{assembly_name}")
    log:
        out = "logs/meryl_{assembly_name}.out",
        err = "logs/meryl_{assembly_name}.err"
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
        db = "results/hifi_reads/meryl/{assembly_name}",
        assembly = "results/{assembly}/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        "results/{assembly}/merqury/{assembly_name}.merqury.qv"
    log:
        out = "logs/merqury_{assembly}_{assembly_name}.out",
        err = "logs/merqury_{assembly}_{assembly_name}.err"
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
        assembly = "results/{assembly}/assembly/{assembly_name}.asm.bp.p_ctg.fa",
        reads = "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
    output:
        directory("results/{assembly}/inspector/{assembly_name}")
    log:
        out = "logs/inspector_{assembly}_{assembly_name}.out",
        err = "logs/inspector_{assembly}_{assembly_name}.err"
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

rule gt_suffixerator:
    input:
        "results/{assembly}/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        assembly = "results/{assembly}/lai/gt_suffixerator/{assembly_name}.fa",
        des = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.des",
        esq = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.esq",
        lcp = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.lcp",
        llv = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.llv",
        md5 = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.md5",
        prj = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.prj",
        sds = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.sds",
        ssp = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.ssp",
        suf = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.suf"
    log:
        out = "logs/gt_suffixerator_{assembly}_{assembly_name}.out",
        err = "logs/gt_suffixerator_{assembly}_{assembly_name}.err"
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
        assembly = "results/{assembly}/lai/gt_suffixerator/{assembly_name}.fa",
        des = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.des",
        esq = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.esq",
        lcp = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.lcp",
        llv = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.llv",
        md5 = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.md5",
        prj = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.prj",
        sds = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.sds",
        ssp = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.ssp",
        suf = "results/{assembly}/lai/gt_suffixerator/{assembly_name}_index.suf"
    output:
        "results/{assembly}/lai/gt_ltrharvest/{assembly_name}.fa.harvest.scn"
    log:
        err = "logs/gt_ltrharvest_{assembly}_{assembly_name}.err"
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
        "results/{assembly}/lai/gt_suffixerator/{assembly_name}.fa"
    output:
        "results/{assembly}/lai/ltr_finder_parallel/{assembly_name}.fa.finder.combine.scn"
    log:
        out = "logs/ltr_finder_parallel_{assembly}_{assembly_name}.out",
        err = "logs/ltr_finder_parallel_{assembly}_{assembly_name}.err"
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
        ltrharvest = "results/{assembly}/lai/gt_ltrharvest/{assembly_name}.fa.harvest.scn",
        ltrfinder = "results/{assembly}/lai/ltr_finder_parallel/{assembly_name}.fa.finder.combine.scn"
    output:
        "results/{assembly}/lai/lib_ltrharvest_ltrfinder/{assembly_name}.fa.rawLTR.scn"
    log:
        err = "logs/merge_ltrharvest_and_ltrfinder_{assembly}_{assembly_name}.err"
    conda:
        "../envs/lai.yml"
    shell:
        "cat {input.ltrharvest} {input.ltrfinder} > {output} 2> {log.err}"

rule ltr_retriever:
    input:
        assembly = "results/{assembly}/lai/gt_suffixerator/{assembly_name}.fa",
        inharvest = "results/{assembly}/lai/lib_ltrharvest_ltrfinder/{assembly_name}.fa.rawLTR.scn"
    output:
        "results/{assembly}/lai/ltr_retriever/{assembly_name}.fa.out.LAI"
    log:
        out = "logs/ltr_retriever_{assembly}_{assembly_name}.out",
        err = "logs/ltr_retriever_{assembly}_{assembly_name}.err"
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
