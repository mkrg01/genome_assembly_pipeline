PMGA_FIGSHARE_ARTICLE_ID = "27201798"
PMGA_FIGSHARE_VERSION = "4"
PMGA_ARCHIVE_NAME = "PMGA.tar.gz"
PMGA_DB = "1"

PGA_V2_COMMIT = "6fa800c2a379ec3356da5ccbf2ea7b903e784f46"
PGA_V2_SCRIPT_SHA256 = "809d04f5070b4d8892fa2707776f47989e0e242ac87c85f6624f5df21c444e9a"
PGA_V2_SCRIPT_URL = (
    "https://raw.githubusercontent.com/quxiaojian/PlastidHub/"
    f"{PGA_V2_COMMIT}/scripts/1.2.PGA_v2.pl"
)


rule download_pmga:
    output:
        archive = f"results/downloads/pmga/v{PMGA_FIGSHARE_VERSION}/{PMGA_ARCHIVE_NAME}",
        bundle = directory(f"results/downloads/pmga/v{PMGA_FIGSHARE_VERSION}/PMGA"),
        manifest = f"results/downloads/pmga/v{PMGA_FIGSHARE_VERSION}/manifest.json"
    log:
        out = "logs/download_pmga.out",
        err = "logs/download_pmga.err"
    conda:
        "../envs/pybase.yml"
    params:
        article_id = PMGA_FIGSHARE_ARTICLE_ID,
        version = PMGA_FIGSHARE_VERSION,
        file_name = PMGA_ARCHIVE_NAME
    shell:
        """
        (
            python3 workflow/scripts/download_pmga.py \
                --article-id {params.article_id:q} \
                --version {params.version:q} \
                --file-name {params.file_name:q} \
                --archive {output.archive:q} \
                --bundle {output.bundle:q} \
                --manifest {output.manifest:q}
        ) > {log.out:q} 2> {log.err:q}
        """


rule download_pga_v2_script:
    output:
        script = f"results/downloads/pga_v2/{PGA_V2_COMMIT}/1.2.PGA_v2.pl",
        manifest = f"results/downloads/pga_v2/{PGA_V2_COMMIT}/manifest.json"
    log:
        out = "logs/download_pga_v2_script.out",
        err = "logs/download_pga_v2_script.err"
    conda:
        "../envs/pybase.yml"
    params:
        url = PGA_V2_SCRIPT_URL,
        sha256 = PGA_V2_SCRIPT_SHA256,
        commit = PGA_V2_COMMIT
    shell:
        """
        (
            python3 workflow/scripts/download_pga_v2_script.py \
                --url {params.url:q} \
                --sha256 {params.sha256:q} \
                --commit {params.commit:q} \
                --output {output.script:q} \
                --manifest {output.manifest:q}
        ) > {log.out:q} 2> {log.err:q}
        """


if "mitochondrion" in configured_oatk_organelles() and configured_organelle_annotation_tool("mitochondrion") == "pmga":
    rule annotate_mitochondrion_pmga:
        input:
            genome = "results/oatk/oatk/{assembly_name}.mito.ctg.fasta",
            pmga_bundle = f"results/downloads/pmga/v{PMGA_FIGSHARE_VERSION}/PMGA"
        output:
            annotation = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/{assembly_name}.mitochondrion.annotation.gbk",
            manifest = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/{assembly_name}.mitochondrion.annotation_manifest.json"
        log:
            out = "logs/annotate_mitochondrion_pmga_{assembly_name}.out",
            err = "logs/annotate_mitochondrion_pmga_{assembly_name}.err"
        conda:
            "../envs/pybase.yml"
        container:
            None
        params:
            db = PMGA_DB
        shell:
            """
            (
                python3 workflow/scripts/run_pmga.py \
                    --pmga-bundle {input.pmga_bundle:q} \
                    --input-fasta {input.genome:q} \
                    --annotation {output.annotation:q} \
                    --manifest {output.manifest:q} \
                    --db {params.db:q} \
                    --prefix {wildcards.assembly_name:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if "mitochondrion" in configured_oatk_organelles() and configured_organelle_annotation_tool("mitochondrion") == "mitoz":
    rule annotate_mitochondrion_mitoz:
        input:
            genome = "results/oatk/oatk/{assembly_name}.mito.ctg.fasta"
        output:
            annotation = "results/organelle_annotation/mitochondrion/mitoz/{assembly_name}/{assembly_name}.mitochondrion.annotation.gbk",
            manifest = "results/organelle_annotation/mitochondrion/mitoz/{assembly_name}/{assembly_name}.mitochondrion.annotation_manifest.json"
        log:
            out = "logs/annotate_mitochondrion_mitoz_{assembly_name}.out",
            err = "logs/annotate_mitochondrion_mitoz_{assembly_name}.err"
        conda:
            "../envs/mitoz.yml"
        threads:
            max(1, int(workflow.cores * 0.5))
        params:
            clade = required_mitoz_config_value("mitoz_clade"),
            genetic_code = required_mitoz_genetic_code()
        shell:
            """
            (
                python3 workflow/scripts/run_mitoz.py \
                    --input-fasta {input.genome:q} \
                    --annotation {output.annotation:q} \
                    --manifest {output.manifest:q} \
                    --prefix {wildcards.assembly_name:q} \
                    --clade {params.clade:q} \
                    --genetic-code {params.genetic_code:q} \
                    --threads {threads}
            ) > {log.out:q} 2> {log.err:q}
            """


if "chloroplast" in configured_oatk_organelles() and configured_organelle_annotation_tool("chloroplast") == "pga_v2":
    rule annotate_chloroplast_pga_v2:
        input:
            genome = "results/oatk/oatk/{assembly_name}.pltd.ctg.fasta",
            script = f"results/downloads/pga_v2/{PGA_V2_COMMIT}/1.2.PGA_v2.pl"
        output:
            annotation = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/{assembly_name}.chloroplast.annotation.gbk",
            manifest = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/{assembly_name}.chloroplast.annotation_manifest.json",
            manual_modification = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/manual_modification.txt"
        log:
            out = "logs/annotate_chloroplast_pga_v2_{assembly_name}.out",
            err = "logs/annotate_chloroplast_pga_v2_{assembly_name}.err"
        conda:
            "../envs/pga_v2.yml"
        params:
            reference_dir = pga_v2_reference_dir,
            form = config.get("pga_v2_form", "circular"),
            ir = config.get("pga_v2_ir", "1000"),
            pidentity = config.get("pga_v2_pidentity", "40"),
            link = config.get("pga_v2_link", "Y"),
            redundancy = config.get("pga_v2_redundancy", "N"),
            qcoverage = config.get("pga_v2_qcoverage", "0.5,2.0"),
            warning = config.get("pga_v2_warning", "warning")
        shell:
            """
            (
                python3 workflow/scripts/run_pga_v2.py \
                    --script {input.script:q} \
                    --input-fasta {input.genome:q} \
                    --reference-dir {params.reference_dir:q} \
                    --annotation {output.annotation:q} \
                    --manifest {output.manifest:q} \
                    --manual-modification {output.manual_modification:q} \
                    --assembly-name {wildcards.assembly_name:q} \
                    --form {params.form:q} \
                    --ir {params.ir:q} \
                    --pidentity {params.pidentity:q} \
                    --link {params.link:q} \
                    --redundancy {params.redundancy:q} \
                    --qcoverage {params.qcoverage:q} \
                    --warning {params.warning:q}
            ) > {log.out:q} 2> {log.err:q}
            """
