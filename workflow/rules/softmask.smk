rule build_repeatmodeler_database:
    input:
        "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
    output:
        nhr = "results/repeatmodeler_db/{sample_id}.nhr",
        nin = "results/repeatmodeler_db/{sample_id}.nin",
        njs = "results/repeatmodeler_db/{sample_id}.njs",
        nnd = "results/repeatmodeler_db/{sample_id}.nnd",
        nni = "results/repeatmodeler_db/{sample_id}.nni",
        nog = "results/repeatmodeler_db/{sample_id}.nog",
        nsq = "results/repeatmodeler_db/{sample_id}.nsq",
        translation = "results/repeatmodeler_db/{sample_id}.translation"
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
        "results/repeatmodeler_db/{sample_id}.nhr"
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
        max(1, int(workflow.cores * 0.9))
    shell:
        """
        (
            RepeatModeler \
                -database $(dirname {input})/$(basename {input} .nhr) \
                -threads {threads} \
                -srand 1 \
                -LTRStruct > {log.out} 2> {log.err}
            rm_dir=$(find . -maxdepth 1 -type d -name "RM_*" | head -n 1)
            if [ -z "$rm_dir" ]; then
                echo "Error: RM_ output directory not found." >&2
                exit 1
            fi
            families_fa=$(find "$rm_dir" -type f -name "{wildcards.sample_id}-families.fa" | head -n 1)
            families_stk=$(find "$rm_dir" -type f -name "{wildcards.sample_id}-families.stk" | head -n 1)
            rmod_log=$(find "$rm_dir" -type f -name "{wildcards.sample_id}-rmod.log" | head -n 1)
            if [ ! -f "$families_fa" ] || [ ! -f "$families_stk" ] || [ ! -f "$rmod_log" ]; then
                echo "Error: Expected output files not found in $rm_dir." >&2
                exit 1
            fi
            mv $rm_dir/* $(dirname {output.consensus})/
            rm -rf $rm_dir
        ) > {log.out} 2> {log.err}
        """

rule download_dfam_database:
    input:
        "results/fcs_gx_db/check"
    output:
        root = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz",
        root_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz.md5",
        lineage = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz",
        lineage_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz.md5"
    log:
        out = "logs/download_dfam_database.out",
        err = "logs/download_dfam_database.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        root = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz",
        root_md5 = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz.md5",
        lineage = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz",
        lineage_md5 = f"https://www.dfam.org/releases/Dfam_{config['dfam_version']}/families/FamDB/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz.md5"
    shell:
        """
        (
            wget -O {output.root} {params.root}
            wget -O {output.root_md5} {params.root_md5}
            wget -O {output.lineage} {params.lineage}
            wget -O {output.lineage_md5} {params.lineage_md5}
            python3 workflow/scripts/dfam_md5_validator.py --file {output.root} --md5_file {output.root_md5}
            python3 workflow/scripts/dfam_md5_validator.py --file {output.lineage} --md5_file {output.lineage_md5}
        ) > {log.out} 2> {log.err}
        """

rule print_dfam_database_info:
    input:
        root = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz",
        root_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz.md5",
        lineage = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz",
        lineage_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz.md5"
    output:
        "results/dfam/dfam_info.txt"
    log:
        "logs/print_dfam_database_info.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "famdb.py \
            -i $(dirname {input.lineage}) \
            info > {output} 2> {log}"

rule print_dfam_repeat_number:
    input:
        root = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz",
        root_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz.md5",
        lineage = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz",
        lineage_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz.md5"
    output:
        f"results/dfam/dfam_repeat_number_{config['dfam_lineage_name']}.txt"
    log:
        "logs/print_dfam_repeat_number.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        lineage_name = config['dfam_lineage_name']
    shell:
        "famdb.py \
            -i $(dirname {input.lineage}) \
            lineage \
            --format totals \
            --ancestors \
            --descendants \
            {params.lineage_name} > {output} 2> {log}"

rule export_dfam_repeat_fasta:
    input:
        root = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz",
        root_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.0.h5.gz.md5",
        lineage = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz",
        lineage_md5 = f"results/dfam/dfam{config['dfam_version'].replace('.', '')}_full.{config['dfam_lineage_id']}.h5.gz.md5"
    output:
        f"results/dfam/dfam_{config['dfam_lineage_name']}.repeat.fasta"
    log:
        "logs/export_dfam_repeat_fasta.err"
    container:
        "docker://dfam/tetools:1.93"
    params:
        lineage_name = config['dfam_lineage_name']
    shell:
        "famdb.py \
            -i $(dirname {input.lineage}) \
            families \
            --format fasta_name \
            --ancestors \
            --descendants \
            --include-class-in-name \
            {params.lineage_name} > {output} 2> {log}"

rule merge_repeat_datasets:
    input:
        repeatmodeler = "results/repeatmodeler/{sample_id}-families.fa",
        dfam = f"results/dfam/dfam_{config['dfam_lineage_name']}.repeat.fasta"
    output:
        "results/repeatmasker_lib/{sample_id}_repeatmasker_lib.fa"
    log:
        "logs/merge_repeat_datasets_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    shell:
        "cat {input.repeatmodeler} {input.dfam} > {output} 2> {log}"

rule repeatmasker:
    input:
        library = "results/repeatmasker_lib/{sample_id}_repeatmasker_lib.fa",
        assembly = "results/fcs_gx_clean/{sample_id}.asm.bp.p_ctg.clean.fa"
    output:
        "results/repeatmasker/{sample_id}.asm.bp.p_ctg.fa.masked"
    log:
        out = "logs/repeatmasker_{sample_id}.out",
        err = "logs/repeatmasker_{sample_id}.err"
    container:
        "docker://dfam/tetools:1.93"
    threads:
        max(1, int(workflow.cores * 0.9))
    shell:
        "RepeatMasker \
            -engine rmblast \
            -parallel $(( {threads} / 4)) \
            -lib {input.library} \
            -dir $(dirname {output}) \
            -xsmall \
            --gff \
            {input.assembly} > {log.out} 2> {log.err}"

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
        python -c 'import sys; num = int(sys.argv[1]); den = int(sys.argv[2]); print("{{:,.1f}}% masked ({{:,}}/{{:,}} bp)".format(num/den*100, num, den))' \
            $num_masked_bp $num_total_bp > {output} 2> {log}
        """
