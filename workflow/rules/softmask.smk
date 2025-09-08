dfam_partitions = config["dfam_partitions"].split(",")

rule build_repeatmodeler_database:
    input:
        "results/fcs/assembly/{sample_id}.asm.bp.p_ctg.fa"
    output:
        nhr = "results/repeatmodeler/{sample_id}.nhr",
        nin = "results/repeatmodeler/{sample_id}.nin",
        njs = "results/repeatmodeler/{sample_id}.njs",
        nnd = "results/repeatmodeler/{sample_id}.nnd",
        nni = "results/repeatmodeler/{sample_id}.nni",
        nog = "results/repeatmodeler/{sample_id}.nog",
        nsq = "results/repeatmodeler/{sample_id}.nsq",
        translation = "results/repeatmodeler/{sample_id}.translation"
    log:
        out = "logs/build_repeatmodeler_database_{sample_id}.out",
        err = "logs/build_repeatmodeler_database_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "BuildDatabase \
            -name $(dirname {output.nhr})/$(basename {output.nhr} .nhr) \
            {input} > {log.out} 2> {log.err}"

rule repeatmodeler:
    input:
        nhr = "results/repeatmodeler/{sample_id}.nhr",
        nin = "results/repeatmodeler/{sample_id}.nin",
        njs = "results/repeatmodeler/{sample_id}.njs",
        nnd = "results/repeatmodeler/{sample_id}.nnd",
        nni = "results/repeatmodeler/{sample_id}.nni",
        nog = "results/repeatmodeler/{sample_id}.nog",
        nsq = "results/repeatmodeler/{sample_id}.nsq",
        translation = "results/repeatmodeler/{sample_id}.translation",
        flag_smudgeplot = "results/hifi_reads/smudgeplot/{sample_id}_masked_errors_smu.txt",
        flag_genomescope = "results/hifi_reads/genomescope2/{sample_id}_linear_plot.png",
        flag_seqkit_hifiasm = "results/hifiasm/seqkit/{sample_id}_seqkit_stats.txt",
        flag_seqkit_fcs = "results/fcs/seqkit/{sample_id}_seqkit_stats.txt",
        flag_busco_hifiasm = "results/hifiasm/busco_genome/BUSCO_{sample_id}.asm.bp.p_ctg.fa",
        flag_busco_fcs = "results/fcs/busco_genome/BUSCO_{sample_id}.asm.bp.p_ctg.fa",
        flag_merqury_hifiasm = "results/hifiasm/merqury/{sample_id}.merqury.qv",
        flag_merqury_fcs = "results/fcs/merqury/{sample_id}.merqury.qv",
        flag_inspector_hifiasm = "results/hifiasm/inspector/{sample_id}",
        flag_inspector_fcs = "results/fcs/inspector/{sample_id}"
    output:
        consensus = "results/repeatmodeler/{sample_id}-families.fa",
        seed = "results/repeatmodeler/{sample_id}-families.stk",
        log = "results/repeatmodeler/{sample_id}-rmod.log"
    log:
        out = "logs/repeatmodeler_{sample_id}.out",
        err = "logs/repeatmodeler_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    threads:
        max(1, int(workflow.cores * 0.8))
    shell:
        """
        (
            cd $(dirname {input.nhr})
            RepeatModeler \
                -database {wildcards.sample_id} \
                -threads {threads} \
                -srand 1 \
                -LTRStruct
            cd ../../
        ) > {log.out} 2> {log.err}
        """

rule download_dfam_database:
    input:
        flag_gxdb = "results/downloads/fcs/.gxdb_checked"
    output:
        db = f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5.gz",
        md5 = f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5.gz.md5"
    log:
        out = "logs/download_dfam_database_{partition}.out",
        err = "logs/download_dfam_database_{partition}.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        url_db = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5.gz",
        url_md5 = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5.gz.md5"
    shell:
        """
        (
            wget -O {output.db} {params.url_db}
            wget -O {output.md5} {params.url_md5}
            python3 workflow/scripts/dfam_md5_validator.py --file {output.db} --md5_file {output.md5}
        ) > {log.out} 2> {log.err}
        """

rule unzip_dfam_database:
    input:
        f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5.gz"
    output:
        f"results/repeatmasker/dfam/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5"
    log:
        out = "logs/unzip_dfam_database_{partition}.out",
        err = "logs/unzip_dfam_database_{partition}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        """
        (
            mkdir -p $(dirname {output})
            zcat {input} > {output}
        ) > {log.out} 2> {log.err}
        """

rule print_dfam_database_info:
    input:
        expand(f"results/repeatmasker/dfam/dfam{config['dfam_version'].replace('.', '')}_full." + {partition} + ".h5", partition=dfam_partitions)
    output:
        "results/repeatmasker/dfam/dfam_info.txt"
    log:
        "logs/print_dfam_database_info.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "famdb.py \
            -i $(dirname {output}) \
            info > {output} 2> {log}"

rule print_dfam_repeat_number:
    input:
        "results/repeatmasker/dfam/dfam_info.txt"
    output:
        f"results/repeatmasker/dfam/dfam_repeat_number_{config['dfam_lineage_name']}.txt"
    log:
        "logs/print_dfam_repeat_number.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        lineage_name = config['dfam_lineage_name']
    shell:
        "famdb.py \
            -i $(dirname {input}) \
            lineage \
            --format totals \
            --ancestors \
            --descendants \
            {params.lineage_name} > {output} 2> {log}"

rule export_dfam_repeat_fasta:
    input:
        f"results/repeatmasker/dfam/dfam_repeat_number_{config['dfam_lineage_name']}.txt"
    output:
        f"results/repeatmasker/dfam/dfam_{config['dfam_lineage_name']}.repeat.fasta"
    log:
        "logs/export_dfam_repeat_fasta.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        lineage_name = config['dfam_lineage_name']
    shell:
        "famdb.py \
            -i $(dirname {input}) \
            families \
            --format fasta_name \
            --ancestors \
            --descendants \
            --include-class-in-name \
            {params.lineage_name} > {output} 2> {log}"

rule merge_repeat_datasets:
    input:
        repeatmodeler = "results/repeatmodeler/{sample_id}-families.fa",
        dfam = f"results/repeatmasker/dfam/dfam_{config['dfam_lineage_name']}.repeat.fasta"
    output:
        "results/repeatmasker/library/{sample_id}_repeatmasker_lib.fa"
    log:
        "logs/merge_repeat_datasets_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "cat {input.repeatmodeler} {input.dfam} > {output} 2> {log}"

rule repeatmasker:
    input:
        library = "results/repeatmasker/library/{sample_id}_repeatmasker_lib.fa",
        assembly = "results/fcs/assembly/{sample_id}.asm.bp.p_ctg.fa",
        flag_lai_1 = "results/hifiasm/lai/ltr_retriever/{sample_id}.fa.out.LAI",
        flag_lai_2 = "results/fcs/lai/ltr_retriever/{sample_id}.fa.out.LAI"
    output:
        "results/repeatmasker/{sample_id}.asm.bp.p_ctg.fa.masked"
    log:
        out = "logs/repeatmasker_{sample_id}.out",
        err = "logs/repeatmasker_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    threads:
        workflow.cores
    shell:
        """
        (
            cd $(dirname {output})
            RepeatMasker \
                -engine rmblast \
                -parallel $(( {threads} / 4)) \
                -lib ../../{input.library} \
                -dir ../../$(dirname {output}) \
                -xsmall \
                --gff \
                ../../{input.assembly}
            cd ../../
            mv .RepeatMaskerCache $(dirname {output})/
        ) > {log.out} 2> {log.err}
        """

rule print_softmasked_percentage:
    input:
        "results/repeatmasker/{sample_id}.asm.bp.p_ctg.fa.masked"
    output:
        "results/repeatmasker/{sample_id}_softmasked_percentage.txt"
    log:
        "logs/print_softmasked_percentage_{sample_id}.log"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        """
        num_total_bp=$(grep -v '^>' {input} | wc -c | tr -d ' ')
        num_masked_bp=$(grep -v '^>' {input} | tr -d -c 'atgc' | wc -c | tr -d ' ')
        python3 -c 'import sys; num = int(sys.argv[1]); den = int(sys.argv[2]); print("{{:,.1f}}% masked ({{:,}}/{{:,}} bp)".format(num/den*100, num, den))' \
            $num_masked_bp $num_total_bp > {output} 2> {log}
        """
