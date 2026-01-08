wildcard_constraints:
    assembly_name = config["assembly_name"]

rule build_repeatmodeler_database:
    input:
        "results/fcs/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        nhr = "results/repeatmodeler/{assembly_name}.nhr",
        nin = "results/repeatmodeler/{assembly_name}.nin",
        njs = "results/repeatmodeler/{assembly_name}.njs",
        nnd = "results/repeatmodeler/{assembly_name}.nnd",
        nni = "results/repeatmodeler/{assembly_name}.nni",
        nog = "results/repeatmodeler/{assembly_name}.nog",
        nsq = "results/repeatmodeler/{assembly_name}.nsq",
        translation = "results/repeatmodeler/{assembly_name}.translation"
    log:
        out = "logs/build_repeatmodeler_database_{assembly_name}.out",
        err = "logs/build_repeatmodeler_database_{assembly_name}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "BuildDatabase \
            -name $(dirname {output.nhr})/$(basename {output.nhr} .nhr) \
            {input} > {log.out} 2> {log.err}"

rule repeatmodeler:
    input:
        nhr = "results/repeatmodeler/{assembly_name}.nhr",
        nin = "results/repeatmodeler/{assembly_name}.nin",
        njs = "results/repeatmodeler/{assembly_name}.njs",
        nnd = "results/repeatmodeler/{assembly_name}.nnd",
        nni = "results/repeatmodeler/{assembly_name}.nni",
        nog = "results/repeatmodeler/{assembly_name}.nog",
        nsq = "results/repeatmodeler/{assembly_name}.nsq",
        translation = "results/repeatmodeler/{assembly_name}.translation"
    output:
        consensus = "results/repeatmodeler/{assembly_name}-families.fa",
        seed = "results/repeatmodeler/{assembly_name}-families.stk",
        log = "results/repeatmodeler/{assembly_name}-rmod.log"
    log:
        out = "logs/repeatmodeler_{assembly_name}.out",
        err = "logs/repeatmodeler_{assembly_name}.err"
    container:
        "docker://dfam/tetools:1.93"
    threads:
        max(1, int(workflow.cores))
    shell:
        """
        (
            cd $(dirname {input.nhr})
            RepeatModeler \
                -database {wildcards.assembly_name} \
                -threads {threads} \
                -srand 1 \
                -LTRStruct
            cd ../../
        ) > {log.out} 2> {log.err}
        """

rule download_dfam_database:
    output:
        db = f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5.gz",
        md5 = f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5.gz.md5"
    log:
        out = "logs/download_dfam_database_{partition}.out",
        err = "logs/download_dfam_database_{partition}.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        url_db = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5.gz",
        url_md5 = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5.gz.md5"
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
        f"results/downloads/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5.gz"
    output:
        f"results/repeatmasker/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5"
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
        expand(f"results/repeatmasker/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{{partition}}.h5", partition=config["dfam_partitions"].split(","))
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

rule export_dfam_repeat_fasta:
    input:
        "results/repeatmasker/dfam/dfam_info.txt"
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
        repeatmodeler = "results/repeatmodeler/{assembly_name}-families.fa",
        dfam = f"results/repeatmasker/dfam/dfam_{config['dfam_lineage_name']}.repeat.fasta"
    output:
        "results/repeatmasker/library/{assembly_name}_repeatmasker_lib.fa"
    log:
        "logs/merge_repeat_datasets_{assembly_name}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "cat {input.repeatmodeler} {input.dfam} > {output} 2> {log}"

rule repeatmasker:
    input:
        library = "results/repeatmasker/library/{assembly_name}_repeatmasker_lib.fa",
        assembly = "results/fcs/assembly/{assembly_name}.asm.bp.p_ctg.fa"
    output:
        masked = "results/repeatmasker/{assembly_name}.asm.bp.p_ctg.fa.masked",
        out_xm = "results/repeatmasker/{assembly_name}.asm.bp.p_ctg.fa.out.xm"
    log:
        out = "logs/repeatmasker_{assembly_name}.out",
        err = "logs/repeatmasker_{assembly_name}.err"
    container:
        "docker://dfam/tetools:1.93"
    threads:
        workflow.cores
    shell:
        """
        (
            cd $(dirname {output.masked})
            RepeatMasker \
                -engine rmblast \
                -parallel $(( {threads} / 4)) \
                -lib ../../{input.library} \
                -dir ../../$(dirname {output.masked}) \
                -xsmall \
                -gff \
                -xm \
                ../../{input.assembly}
            cd ../../
            mv .RepeatMaskerCache $(dirname {output.masked})/
        ) > {log.out} 2> {log.err}
        """
