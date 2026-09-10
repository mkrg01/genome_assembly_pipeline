TE_TOOLS_IMAGE = "docker://dfam/tetools:2.00"
DFAM_COMPONENT_PARTITIONS = (
    "0",
    "curated.consensus.0",
    "uncurated.consensus.0",
    "uncurated.consensus.1",
)

dfam_version = config["dfam_version"]
if dfam_version != "4.0":
    raise ValueError('dfam_version must be "4.0" for the FamDB 3 consensus layout.')
dfam_lineage = config["dfam_lineage_name"]
if not isinstance(dfam_lineage, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9 ._-]*", dfam_lineage):
    raise ValueError("dfam_lineage_name must be a taxon name or taxonomy ID.")

dfam_prefix = f"dfam{dfam_version.replace('.', '')}"
dfam_download_dir = f"results/downloads/dfam/{dfam_version}"
dfam_database_dir = f"results/repeatmasker/dfam/{dfam_version}"
dfam_database_files = expand(
    f"{dfam_database_dir}/famdb/{dfam_prefix}.{{component_partition}}.h5",
    component_partition=DFAM_COMPONENT_PARTITIONS,
)
dfam_repeat_fasta = f"{dfam_database_dir}/dfam_{dfam_lineage}.repeat.fasta"

wildcard_constraints:
    assembly_name = organism_name,
    selected_assembly = selected_assembly_pattern,
    component_partition = "(?:" + "|".join(re.escape(part) for part in DFAM_COMPONENT_PARTITIONS) + ")"

rule build_repeatmodeler_database:
    input:
        lambda wildcards: downstream_assembly_path(
            wildcards.assembly_name,
            wildcards.selected_assembly,
        )
    output:
        database = directory("results/repeatmodeler/database/{selected_assembly}/{assembly_name}")
    log:
        out = "logs/build_repeatmodeler_database_{selected_assembly}_{assembly_name}.out",
        err = "logs/build_repeatmodeler_database_{selected_assembly}_{assembly_name}.err"
    container:
        TE_TOOLS_IMAGE
    shell:
        """
        (
            mkdir -p {output.database:q}
            BuildDatabase -name {output.database:q}/{wildcards.assembly_name:q} {input:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule repeatmodeler:
    input:
        database = "results/repeatmodeler/database/{selected_assembly}/{assembly_name}"
    output:
        consensus = "results/repeatmodeler/{selected_assembly}/{assembly_name}-families.fa",
        seed = "results/repeatmodeler/{selected_assembly}/{assembly_name}-families.stk",
        log = "results/repeatmodeler/{selected_assembly}/{assembly_name}-rmod.log"
    log:
        out = "logs/repeatmodeler_{selected_assembly}_{assembly_name}.out",
        err = "logs/repeatmodeler_{selected_assembly}_{assembly_name}.err"
    container:
        TE_TOOLS_IMAGE
    threads:
        max(1, int(workflow.cores))
    shell:
        """
        (
            database=$(realpath {input.database:q})/{wildcards.assembly_name:q}
            mkdir -p "$(dirname {output.consensus:q})"
            consensus=$(realpath -m {output.consensus:q})
            seed=$(realpath -m {output.seed:q})
            modeler_log=$(realpath -m {output.log:q})
            rm -f -- "$database-families.fa" "$database-families.stk" "$database-rmod.log"
            cd "$(dirname "$consensus")"
            RepeatModeler \
                -database "$database" \
                -threads {threads} \
                -srand 1 \
                -LTRStruct
            mv -- "$database-families.fa" "$consensus"
            mv -- "$database-families.stk" "$seed"
            mv -- "$database-rmod.log" "$modeler_log"
        ) > {log.out:q} 2> {log.err:q}
        """

rule download_dfam_database:
    input:
        validator = os.path.join(workflow.basedir, "scripts/dfam_md5_validator.py")
    output:
        db = f"{dfam_download_dir}/{dfam_prefix}.{{component_partition}}.h5.gz",
        md5 = f"{dfam_download_dir}/{dfam_prefix}.{{component_partition}}.h5.gz.md5"
    log:
        out = "logs/download_dfam_database_{component_partition}.out",
        err = "logs/download_dfam_database_{component_partition}.err"
    container:
        TE_TOOLS_IMAGE
    params:
        url = f"https://www.dfam.org/releases/Dfam_{dfam_version}/families/FamDB/{dfam_prefix}.{{component_partition}}.h5.gz"
    shell:
        """
        (
            wget -O {output.db:q} {params.url:q}
            wget -O {output.md5:q} {params.url:q}.md5
            python3 {input.validator:q} --file {output.db:q} --md5_file {output.md5:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule unzip_dfam_database:
    input:
        db = f"{dfam_download_dir}/{dfam_prefix}.{{component_partition}}.h5.gz",
        md5 = f"{dfam_download_dir}/{dfam_prefix}.{{component_partition}}.h5.gz.md5"
    output:
        f"{dfam_database_dir}/famdb/{dfam_prefix}.{{component_partition}}.h5"
    log:
        out = "logs/unzip_dfam_database_{component_partition}.out",
        err = "logs/unzip_dfam_database_{component_partition}.err"
    container:
        TE_TOOLS_IMAGE
    shell:
        """
        (
            mkdir -p "$(dirname {output:q})"
            zcat {input.db:q} > {output:q}
        ) > {log.out:q} 2> {log.err:q}
        """

rule print_dfam_database_info:
    input:
        dfam_database_files
    output:
        f"{dfam_database_dir}/dfam_info.txt"
    log:
        "logs/print_dfam_database_info.err"
    container:
        TE_TOOLS_IMAGE
    params:
        database = lambda wildcards, input: os.path.dirname(input[0])
    shell:
        "python3 /opt/FamDB/famdb.py \
            -i {params.database:q} \
            info > {output:q} 2> {log:q}"

rule export_dfam_repeat_fasta:
    input:
        database = dfam_database_files,
        info = f"{dfam_database_dir}/dfam_info.txt"
    output:
        dfam_repeat_fasta
    log:
        "logs/export_dfam_repeat_fasta.err"
    container:
        TE_TOOLS_IMAGE
    params:
        database = lambda wildcards, input: os.path.dirname(input.database[0]),
        lineage_name = dfam_lineage
    shell:
        """
        (
            python3 /opt/FamDB/famdb.py \
                -i {params.database:q} \
                families \
                --format fasta_name \
                --ancestors \
                --descendants \
                --include-class-in-name \
                {params.lineage_name:q} > {output:q}
            if ! grep -q '^>' {output:q}; then
                echo "No Dfam consensus sequences exported; check dfam_lineage_name." >&2
                exit 1
            fi
        ) 2> {log:q}
        """

rule merge_repeat_datasets:
    input:
        repeatmodeler = "results/repeatmodeler/{selected_assembly}/{assembly_name}-families.fa",
        dfam = dfam_repeat_fasta
    output:
        "results/repeatmasker/library/{selected_assembly}/{assembly_name}_repeatmasker_lib.fa"
    log:
        "logs/merge_repeat_datasets_{selected_assembly}_{assembly_name}.err"
    container:
        TE_TOOLS_IMAGE
    shell:
        "cat {input.repeatmodeler:q} {input.dfam:q} > {output:q} 2> {log:q}"

rule repeatmasker:
    input:
        library = "results/repeatmasker/library/{selected_assembly}/{assembly_name}_repeatmasker_lib.fa",
        assembly = lambda wildcards: downstream_assembly_path(
            wildcards.assembly_name,
            wildcards.selected_assembly,
        )
    output:
        masked = "results/repeatmasker/{selected_assembly}/{assembly_name}.fa.masked",
        out_xm = "results/repeatmasker/{selected_assembly}/{assembly_name}.fa.out.xm"
    log:
        out = "logs/repeatmasker_{selected_assembly}_{assembly_name}.out",
        err = "logs/repeatmasker_{selected_assembly}_{assembly_name}.err"
    container:
        TE_TOOLS_IMAGE
    threads:
        workflow.cores
    shell:
        """
        (
            library=$(realpath {input.library:q})
            assembly=$(realpath {input.assembly:q})
            output_dir=$(dirname {output.masked:q})
            mkdir -p "$output_dir"
            cd "$output_dir"
            batches=$(( {threads} / 4 ))
            if [ "$batches" -lt 1 ]; then batches=1; fi
            RepeatMasker \
                -engine rmblast \
                -parallel "$batches" \
                -lib "$library" \
                -dir "$PWD" \
                -xsmall \
                -gff \
                -xm \
                "$assembly"
        ) > {log.out:q} 2> {log.err:q}
        """
