import os
import glob
rnaseq_samples = glob.glob("raw_data/*_1.fastq.gz")
if not rnaseq_samples:
    raise RuntimeError("No RNA-Seq files matching 'raw_data/*_1.fastq.gz' were found. Please check the input directory.")
rnaseq_sample_ids = sorted({os.path.basename(f).replace("_1.fastq.gz", "") for f in rnaseq_samples})

wildcard_constraints:
    assembly_name = config["assembly_name"]

rule fastp_rnaseq_sample:
    input:
        rnaseq_1 = "raw_data/{rnaseq_sample_id}_1.fastq.gz",
        rnaseq_2 = "raw_data/{rnaseq_sample_id}_2.fastq.gz"
    output:
        rnaseq_1 = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_1.fastq",
        rnaseq_2 = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_2.fastq",
        html = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_fastp.html",
        json = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_fastp.json"
    log:
        out = "logs/fastp_rnaseq_sample_{rnaseq_sample_id}.out",
        err = "logs/fastp_rnaseq_sample_{rnaseq_sample_id}.err"
    conda:
        "../envs/fastp.yml"
    threads:
        max(1, int(workflow.cores * 0.1))
    shell:
        "fastp \
            --in1 {input.rnaseq_1} \
            --in2 {input.rnaseq_2} \
            --out1 {output.rnaseq_1} \
            --out2 {output.rnaseq_2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} > {log.out} 2> {log.err}"

rule download_orthodb_proteins:
    output:
        f"results/downloads/orthodb/{config['orthodb_lineage']}.fa"
    log:
        out = "logs/download_orthodb_proteins.out",
        err = "logs/download_orthodb_proteins.err"
    container:
        "docker://teambraker/braker3:v3.0.7.6"
    params:
        url = f"https://bioinf.uni-greifswald.de/bioinf/partitioned_odb{config['orthodb_version']}/{config['orthodb_lineage']}.fa.gz",
        md5sum = config['orthodb_md5sum']
    shell:
        """
        (
            wget -O {output}.gz {params.url}
            actual_md5=$(md5sum {output}.gz | cut -d' ' -f1)
            if [ "$actual_md5" != "{params.md5sum}" ]; then
                echo "MD5 checksum mismatch for {output}.gz"
                echo "Expected: {params.md5sum}, but got: $actual_md5"
                echo "Please check the download URL or the MD5 checksum."
                echo "Exiting with error code 1."
                exit 1
            else
                echo "MD5 checksum verified for {output}.gz"
            fi
            gunzip {output}.gz
        ) > {log.out} 2> {log.err}
        """

rule braker3:
    input:
        assembly = "results/repeatmasker/{assembly_name}.asm.bp.p_ctg.fa.masked",
        rnaseq = expand("results/rnaseq_reads/fastp/{rnaseq_sample_id}_{pair}.fastq", rnaseq_sample_id=rnaseq_sample_ids, pair=[1, 2]),
        protein_dataset = f"results/downloads/orthodb/{config['orthodb_lineage']}.fa"
    output:
        gff3 = "results/braker3/{assembly_name}/braker.gff3",
        gtf = "results/braker3/{assembly_name}/braker.gtf",
        cds = "results/braker3/{assembly_name}/braker.codingseq",
        aa = "results/braker3/{assembly_name}/braker.aa",
        augustus_config = directory("results/braker3/{assembly_name}/augustus_config")
    log:
        out = "logs/braker3_{assembly_name}.out",
        err = "logs/braker3_{assembly_name}.err"
    container:
        "docker://teambraker/braker3:v3.0.7.6"
    threads:
        48
    params:
        rnaseq_ids=",".join(rnaseq_sample_ids),
        rnaseq_dir=lambda wildcards, input: os.path.commonpath(input.rnaseq)
    shell:
        """
        (
            cp -r /opt/Augustus/config {output.augustus_config} # https://github.com/Gaius-Augustus/BRAKER/issues/609
            export AUGUSTUS_CONFIG_PATH=$PWD/{output.augustus_config}
            braker.pl \
                --species={wildcards.assembly_name} \
                --genome={input.assembly} \
                --prot_seq={input.protein_dataset} \
                --rnaseq_sets_ids={params.rnaseq_ids} \
                --rnaseq_sets_dirs={params.rnaseq_dir} \
                --workingdir=$(dirname {output.gff3}) \
                --gff3 \
                --threads {threads}
        ) > {log.out} 2> {log.err}
        """

rule copy_isoforms:
    input:
        cds = "results/braker3/{assembly_name}/braker.codingseq",
        aa = "results/braker3/{assembly_name}/braker.aa"
    output:
        cds = "results/isoforms/{assembly_name}_cds.fa",
        aa = "results/isoforms/{assembly_name}_aa.fa"
    log:
        out = "logs/copy_isoforms_{assembly_name}.out",
        err = "logs/copy_isoforms_{assembly_name}.err"
    conda:
        "../envs/cdskit.yml"
    shell:
        """
        (
            cp {input.cds} {output.cds}
            cp {input.aa} {output.aa}
        ) > {log.out} 2> {log.err}
        """

rule extract_longest_cds:
    input:
        cds = "results/braker3/{assembly_name}/braker.codingseq",
        aa = "results/braker3/{assembly_name}/braker.aa"
    output:
        cds = "results/longest_cds/{assembly_name}_cds.fa",
        aa = "results/longest_cds/{assembly_name}_aa.fa"
    log:
        out = "logs/extract_longest_cds_{assembly_name}.out",
        err = "logs/extract_longest_cds_{assembly_name}.err"
    conda:
        "../envs/cdskit.yml"
    threads:
        4
    shell:
        """
        (
            seqkit seq --threads {threads} {input.cds} \
            | cdskit aggregate --expression "t.*" \
            | seqkit seq --threads {threads} --out-file {output.cds}

            seqkit seq --threads {threads} {input.aa} \
            | cdskit aggregate --expression "t.*" \
            | seqkit seq --threads {threads} --out-file {output.aa}
        ) > {log.out} 2> {log.err}
        """

rule seqkit_stats_proteins:
    input:
        "results/{gene}/{assembly_name}_aa.fa"
    output:
        txt = "results/{gene}/seqkit/{assembly_name}_seqkit_stats.txt",
        tsv = "results/{gene}/seqkit/{assembly_name}_seqkit_stats.tsv"
    log:
        out = "logs/seqkit_stats_proteins_{gene}_{assembly_name}.out",
        err = "logs/seqkit_stats_proteins_{gene}_{assembly_name}.err"
    conda:
        "../envs/seqkit.yml"
    wildcard_constraints:
        gene = "(isoforms|longest_cds)"
    shell:
        """
        (
            seqkit stats --all {input} > {output.txt}
            seqkit stats --all --tabular {input} > {output.tsv}
        ) > {log.out} 2> {log.err}
        """

rule busco_proteins_mode:
    input:
        aa = "results/{gene}/{assembly_name}_aa.fa",
        database = "results/downloads/busco_downloads"
    output:
        directory("results/{gene}/busco_proteins/BUSCO_{assembly_name}_aa.fa")
    log:
        out = "logs/busco_proteins_mode_{gene}_{assembly_name}.out",
        err = "logs/busco_proteins_mode_{gene}_{assembly_name}.err"
    conda:
        "../envs/busco.yml"
    threads:
        max(1, int(workflow.cores * 0.95))
    params:
        lineage_dataset = config["busco_lineage_dataset"]
    wildcard_constraints:
        gene = "(isoforms|longest_cds)"
    shell:
        "busco \
            --in {input.aa} \
            --out_path $(dirname {output}) \
            --mode proteins \
            --cpu {threads} \
            --lineage_dataset {params.lineage_dataset} \
            --download_path {input.database} \
            --offline > {log.out} 2> {log.err}"

rule download_omamer_database:
    output:
        "results/downloads/omamer/LUCA.h5"
    log:
        out = "logs/download_omamer_database.out",
        err = "logs/download_omamer_database.err"
    conda:
        "../envs/omamer.yml"
    params:
        url = "https://omabrowser.org/All/LUCA.h5"
    shell:
        "wget {params.url} -O {output} > {log.out} 2> {log.err}"

rule omamer_search:
    input:
        aa = "results/longest_cds/{assembly_name}_aa.fa",
        database = "results/downloads/omamer/LUCA.h5"
    output:
        "results/longest_cds/omark/{assembly_name}.omamer"
    log:
        out = "logs/omamer_search_{assembly_name}.out",
        err = "logs/omamer_search_{assembly_name}.err"
    conda:
        "../envs/omamer.yml"
    threads:
        max(1, int(workflow.cores * 0.95))
    shell:
        "omamer search \
            --db {input.database} \
            --query {input.aa} \
            --out {output} \
            --nthreads {threads} > {log.out} 2> {log.err}"

rule omark:
    input:
        database = "results/downloads/omamer/LUCA.h5",
        omamer_search = "results/longest_cds/omark/{assembly_name}.omamer"
    output:
        directory("results/longest_cds/omark/{assembly_name}_omark")
    log:
        out = "logs/omark_{assembly_name}.out",
        err = "logs/omark_{assembly_name}.err"
    conda:
        "../envs/omark.yml"
    shell:
        "omark \
            --file {input.omamer_search} \
            --database {input.database} \
            --outputFolder {output} > {log.out} 2> {log.err}"
