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


if configured_oatk_organelles():
    rule prefix_organelle_fasta_record_ids:
        input:
            lambda wildcards: organelle_oatk_genome_path(
                wildcards.assembly_name,
                wildcards.organelle,
            )
        output:
            fasta = (
                "results/organelle_annotation/{organelle}/prefixed/{assembly_name}/"
                "{assembly_name}.{organelle}.ctg.fasta"
            ),
            manifest = (
                "results/organelle_annotation/{organelle}/prefixed/{assembly_name}/"
                "{assembly_name}.{organelle}.ctg.manifest.json"
            )
        log:
            out = "logs/prefix_organelle_fasta_record_ids_{organelle}_{assembly_name}.out",
            err = "logs/prefix_organelle_fasta_record_ids_{organelle}_{assembly_name}.err"
        conda:
            "../envs/pybase.yml"
        wildcard_constraints:
            organelle = "|".join(configured_oatk_organelles())
        params:
            prefix = lambda wildcards: ORGANELLE_FASTA_ID_PREFIXES[
                normalize_organelle_name("oatk_organelle", wildcards.organelle)
            ]
        shell:
            """
            (
                python3 workflow/scripts/prefix_fasta_record_ids.py \
                    --input {input:q} \
                    --output {output.fasta:q} \
                    --manifest {output.manifest:q} \
                    --prefix {params.prefix:q} \
                    --label {wildcards.organelle:q}
            ) > {log.out:q} 2> {log.err:q}
            """


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
            genome = organelle_prefixed_genome_path("{assembly_name}", "mitochondrion"),
            pmga_bundle = f"results/downloads/pmga/v{PMGA_FIGSHARE_VERSION}/PMGA"
        output:
            annotation = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/{assembly_name}.mitochondrion.pre_rna_editing.gbk",
            genome = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/{assembly_name}.mitochondrion.ctg.annotation.fasta",
            manifest = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/{assembly_name}.mitochondrion.manifest.json",
            post_curation = "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/post_curation.pre_rna_editing.md"
        log:
            out = "logs/annotate_mitochondrion_pmga_{assembly_name}.out",
            err = "logs/annotate_mitochondrion_pmga_{assembly_name}.err"
        conda:
            "../envs/pybase.yml"
        container:
            None
        params:
            db = PMGA_DB,
            taxid = taxid or ""
        shell:
            """
            (
                python3 workflow/scripts/run_pmga.py \
                    --pmga-bundle {input.pmga_bundle:q} \
                    --input-fasta {input.genome:q} \
                    --annotation {output.annotation:q} \
                    --annotation-fasta {output.genome:q} \
                    --manifest {output.manifest:q} \
                    --post-curation {output.post_curation:q} \
                    --db {params.db:q} \
                    --prefix {wildcards.assembly_name:q} \
                    --assembly-name {wildcards.assembly_name:q} \
                    --taxid {params.taxid:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if "chloroplast" in configured_oatk_organelles() and configured_organelle_annotation_tool("chloroplast") == "pga_v2":
    rule annotate_chloroplast_pga_v2:
        input:
            genome = organelle_prefixed_genome_path("{assembly_name}", "chloroplast"),
            script = f"results/downloads/pga_v2/{PGA_V2_COMMIT}/1.2.PGA_v2.pl",
            hifi_reads = pga_v2_hifi_reads_for_sequence_fix
        output:
            annotation = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/{assembly_name}.chloroplast.pre_rna_editing.gbk",
            genome = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/{assembly_name}.chloroplast.ctg.annotation.fasta",
            manifest = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/{assembly_name}.chloroplast.manifest.json",
            post_curation = "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/post_curation.pre_rna_editing.md",
            reference_cds_qc_pre = (
                "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/"
                "{assembly_name}.chloroplast.reference_cds_qc.pre_rna_editing.tsv"
            ),
            reference_cds_frameshift_candidates = (
                "results/organelle_annotation/chloroplast/pga_v2/{assembly_name}/"
                "{assembly_name}.chloroplast.reference_cds_qc.frameshift_candidates.json"
            )
        log:
            out = "logs/annotate_chloroplast_pga_v2_{assembly_name}.out",
            err = "logs/annotate_chloroplast_pga_v2_{assembly_name}.err"
        conda:
            "../envs/pga_v2.yml"
        threads:
            max(1, int(workflow.cores * 0.5))
        params:
            reference_dir = required_chloroplast_reference_cds_qc_dir(),
            form = config.get("pga_v2_form", "circular"),
            ir = config.get("pga_v2_ir", "1000"),
            pidentity = config.get("pga_v2_pidentity", "40"),
            link = config.get("pga_v2_link", "Y"),
            redundancy = config.get("pga_v2_redundancy", "N"),
            qcoverage = config.get("pga_v2_qcoverage", "0.5,2.0"),
            warning = config.get("pga_v2_warning", "warning"),
            taxid = taxid or "",
            hifi_reads_arg = pga_v2_hifi_reads_arg,
            fix_hifi_frameshifts_arg = (
                "--fix-hifi-frameshifts"
                if organelle_reference_cds_qc_fix_hifi_frameshifts("chloroplast")
                else ""
            )
        shell:
            """
            (
                python3 workflow/scripts/run_pga_v2.py \
                    --script {input.script:q} \
                    --input-fasta {input.genome:q} \
                    --reference-dir {params.reference_dir:q} \
                    --annotation {output.annotation:q} \
                    --annotation-fasta {output.genome:q} \
                    --manifest {output.manifest:q} \
                    --post-curation {output.post_curation:q} \
                    --reference-cds-qc-pre {output.reference_cds_qc_pre:q} \
                    --reference-cds-frameshift-candidates {output.reference_cds_frameshift_candidates:q} \
                    --assembly-name {wildcards.assembly_name:q} \
                    --form {params.form:q} \
                    --ir {params.ir:q} \
                    --pidentity {params.pidentity:q} \
                    --link {params.link:q} \
                    --redundancy {params.redundancy:q} \
                    --qcoverage {params.qcoverage:q} \
                    --warning {params.warning:q} \
                    --taxid {params.taxid:q} \
                    --threads {threads} \
                    {params.hifi_reads_arg} \
                    {params.fix_hifi_frameshifts_arg}
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_rna_editing_post_curation():
    rule build_organelle_rna_editing_reference:
        input:
            unpack(
                lambda wildcards: organelle_rna_editing_reference_input_paths(
                    wildcards.assembly_name
                )
            )
        output:
            **organelle_rna_editing_reference_paths("{assembly_name}")
        log:
            out = "logs/build_organelle_rna_editing_reference_{assembly_name}.out",
            err = "logs/build_organelle_rna_editing_reference_{assembly_name}.err"
        conda:
            "../envs/organelle_rna_editing.yml"
        params:
            organelles = lambda wildcards, input: " ".join(
                f"--organelle {name}={shlex.quote(str(input[key]))}"
                for key, name in (
                    ("mito", "mitochondrion"),
                    ("pltd", "chloroplast"),
                )
                if key in input.keys()
            )
        shell:
            """
            (
                python3 workflow/scripts/build_organelle_rna_editing_reference.py \
                    --nuclear {input.nuclear:q} \
                    {params.organelles} \
                    --reference {output.reference:q} \
                    --manifest {output.manifest:q}
                samtools faidx {output.reference:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_rna_editing_post_curation():
    rule map_hifi_reads_to_organelle_rna_editing_reference:
        input:
            reference = organelle_rna_editing_reference_paths("{assembly_name}")["reference"],
            reads = "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz"
        output:
            **organelle_rna_editing_hifi_bam_paths("{assembly_name}")
        log:
            out = "logs/map_hifi_reads_to_organelle_rna_editing_reference_{assembly_name}.out",
            err = "logs/map_hifi_reads_to_organelle_rna_editing_reference_{assembly_name}.err"
        conda:
            "../envs/organelle_rna_editing.yml"
        threads:
            max(1, int(workflow.cores * 0.5))
        shell:
            """
            (
                minimap2 \
                    -a \
                    -x map-hifi \
                    {input.reference:q} \
                    {input.reads:q} \
                    -t {threads} \
                | samtools sort \
                    -m 2G \
                    -@ {threads} \
                    -o {output.bam:q}
                samtools index {output.bam:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_rna_editing_post_curation():
    rule map_rnaseq_reads_to_organelle_rna_editing_reference:
        input:
            reference = organelle_rna_editing_reference_paths("{assembly_name}")["reference"],
            rnaseq_1 = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_1.fastq",
            rnaseq_2 = "results/rnaseq_reads/fastp/{rnaseq_sample_id}_2.fastq"
        output:
            bam = (
                "results/organelle_annotation/rna_editing/{assembly_name}/rnaseq/"
                "{rnaseq_sample_id}.rnaseq_to_nuclear_organelle.bam"
            ),
            bai = (
                "results/organelle_annotation/rna_editing/{assembly_name}/rnaseq/"
                "{rnaseq_sample_id}.rnaseq_to_nuclear_organelle.bam.bai"
            )
        log:
            out = (
                "logs/map_rnaseq_reads_to_organelle_rna_editing_reference_"
                "{rnaseq_sample_id}_{assembly_name}.out"
            ),
            err = (
                "logs/map_rnaseq_reads_to_organelle_rna_editing_reference_"
                "{rnaseq_sample_id}_{assembly_name}.err"
            )
        conda:
            "../envs/organelle_rna_editing.yml"
        threads:
            max(1, int(workflow.cores * 0.3))
        wildcard_constraints:
            rnaseq_sample_id = r".+"
        shell:
            """
            (
                minimap2 \
                    -a \
                    -x sr \
                    {input.reference:q} \
                    {input.rnaseq_1:q} \
                    {input.rnaseq_2:q} \
                    -t {threads} \
                | samtools sort \
                    -m 2G \
                    -@ {threads} \
                    -o {output.bam:q}
                samtools index {output.bam:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_rna_editing_post_curation():
    rule curate_organelle_rna_editing:
        input:
            annotation = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.pre_rna_editing.gbk"
            ),
            post_curation = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "post_curation.pre_rna_editing.md"
            ),
            rna_bams = organelle_rna_editing_rnaseq_bam_paths,
            rna_bais = lambda wildcards: [
                f"{path}.bai"
                for path in organelle_rna_editing_rnaseq_bam_paths(wildcards)
            ],
            dna_bam = organelle_rna_editing_hifi_bam_paths("{assembly_name}")["bam"],
            dna_bai = organelle_rna_editing_hifi_bam_paths("{assembly_name}")["bai"]
        output:
            annotation = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.gbk"
            ),
            evidence = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.rna_editing.evidence.tsv"
            ),
            decisions = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.rna_editing.decisions.json"
            ),
            summary = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.rna_editing.summary.tsv"
            ),
            post_curation = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "post_curation.md"
            )
        log:
            out = "logs/curate_organelle_rna_editing_{organelle}_{tool}_{assembly_name}.out",
            err = "logs/curate_organelle_rna_editing_{organelle}_{tool}_{assembly_name}.err"
        conda:
            "../envs/organelle_rna_editing.yml"
        wildcard_constraints:
            organelle = "mitochondrion|chloroplast",
            tool = "pmga|pga_v2"
        params:
            transl_table = lambda wildcards: "11" if wildcards.organelle == "chloroplast" else "1",
            min_rna_depth = ORGANELLE_RNA_EDITING_THRESHOLDS["min_rna_depth"],
            min_edited_reads = ORGANELLE_RNA_EDITING_THRESHOLDS["min_edited_reads"],
            min_edit_fraction = ORGANELLE_RNA_EDITING_THRESHOLDS["min_edit_fraction"],
            essential_rescue_min_rna_depth = ORGANELLE_RNA_EDITING_THRESHOLDS["essential_rescue_min_rna_depth"],
            essential_rescue_min_edited_reads = ORGANELLE_RNA_EDITING_THRESHOLDS["essential_rescue_min_edited_reads"],
            essential_rescue_min_edit_fraction = ORGANELLE_RNA_EDITING_THRESHOLDS["essential_rescue_min_edit_fraction"],
            essential_rescue_max_dna_alt_fraction = ORGANELLE_RNA_EDITING_THRESHOLDS["essential_rescue_max_dna_alt_fraction"],
            min_base_quality = ORGANELLE_RNA_EDITING_THRESHOLDS["min_base_quality"],
            min_mapping_quality = ORGANELLE_RNA_EDITING_THRESHOLDS["min_mapping_quality"],
            min_dna_depth = ORGANELLE_RNA_EDITING_THRESHOLDS["min_dna_depth"],
            max_dna_alt_fraction = ORGANELLE_RNA_EDITING_THRESHOLDS["max_dna_alt_fraction"],
            reference_dir_arg = lambda wildcards: organelle_reference_cds_qc_reference_dir_arg(wildcards.organelle)
        shell:
            """
            (
                python3 workflow/scripts/curate_organelle_rna_editing.py \
                    --input-gbk {input.annotation:q} \
                    --output-gbk {output.annotation:q} \
                    --evidence-tsv {output.evidence:q} \
                    --decisions-json {output.decisions:q} \
                    --summary-tsv {output.summary:q} \
                    --input-post-curation {input.post_curation:q} \
                    --output-post-curation {output.post_curation:q} \
                    --organelle {wildcards.organelle:q} \
                    --tool {wildcards.tool:q} \
                    --transl-table {params.transl_table:q} \
                    --rna-bam {input.rna_bams:q} \
                    --dna-bam {input.dna_bam:q} \
                    --min-rna-depth {params.min_rna_depth} \
                    --min-edited-reads {params.min_edited_reads} \
                    --min-edit-fraction {params.min_edit_fraction} \
                    --essential-rescue-min-rna-depth {params.essential_rescue_min_rna_depth} \
                    --essential-rescue-min-edited-reads {params.essential_rescue_min_edited_reads} \
                    --essential-rescue-min-edit-fraction {params.essential_rescue_min_edit_fraction} \
                    --essential-rescue-max-dna-alt-fraction {params.essential_rescue_max_dna_alt_fraction} \
                    --min-base-quality {params.min_base_quality} \
                    --min-mapping-quality {params.min_mapping_quality} \
                    --min-dna-depth {params.min_dna_depth} \
                    --max-dna-alt-fraction {params.max_dna_alt_fraction} \
                    {params.reference_dir_arg}
            ) > {log.out:q} 2> {log.err:q}
            """


if "mitochondrion" in configured_oatk_organelles_with_annotation() and configured_organelle_annotation_tool("mitochondrion") == "pmga":
    rule reference_cds_qc_pre_rna_editing_pmga:
        input:
            annotation = (
                "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/"
                "{assembly_name}.mitochondrion.pre_rna_editing.gbk"
            )
        output:
            qc_tsv = (
                "results/organelle_annotation/mitochondrion/pmga/{assembly_name}/"
                "{assembly_name}.mitochondrion.reference_cds_qc.pre_rna_editing.tsv"
            )
        log:
            out = "logs/reference_cds_qc_pre_rna_editing_mitochondrion_pmga_{assembly_name}.out",
            err = "logs/reference_cds_qc_pre_rna_editing_mitochondrion_pmga_{assembly_name}.err"
        conda:
            "../envs/pybase.yml"
        params:
            reference_dir = lambda wildcards: organelle_reference_cds_qc_reference_dir("mitochondrion")
        shell:
            """
            (
                python3 workflow/scripts/reference_cds_qc.py \
                    --annotation {input.annotation:q} \
                    --reference-dir {params.reference_dir:q} \
                    --qc-tsv {output.qc_tsv:q} \
                    --organelle mitochondrion \
                    --tool pmga \
                    --phase pre_rna_editing \
                    --default-transl-table 1
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_rna_editing_post_curation():
    rule reference_cds_qc_post_rna_editing:
        input:
            annotation = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.gbk"
            )
        output:
            qc_tsv = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.reference_cds_qc.post_rna_editing.tsv"
            ),
            manual_candidates = (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.manual_rna_editing_candidates.tsv"
            )
        log:
            out = "logs/reference_cds_qc_post_rna_editing_{organelle}_{tool}_{assembly_name}.out",
            err = "logs/reference_cds_qc_post_rna_editing_{organelle}_{tool}_{assembly_name}.err"
        conda:
            "../envs/pybase.yml"
        wildcard_constraints:
            organelle = "mitochondrion|chloroplast",
            tool = "pmga|pga_v2"
        params:
            reference_dir = lambda wildcards: organelle_reference_cds_qc_reference_dir(wildcards.organelle),
            transl_table = lambda wildcards: "11" if wildcards.organelle == "chloroplast" else "1"
        shell:
            """
            (
                python3 workflow/scripts/reference_cds_qc.py \
                    --annotation {input.annotation:q} \
                    --reference-dir {params.reference_dir:q} \
                    --qc-tsv {output.qc_tsv:q} \
                    --manual-candidates-tsv {output.manual_candidates:q} \
                    --organelle {wildcards.organelle:q} \
                    --tool {wildcards.tool:q} \
                    --phase post_rna_editing \
                    --default-transl-table {params.transl_table:q}
            ) > {log.out:q} 2> {log.err:q}
            """


if configured_oatk_organelles_with_annotation():
    rule plot_organelle_pycirclize:
        input:
            annotation = (
                lambda wildcards: organelle_annotation_input_for_downstream(
                    wildcards.assembly_name,
                    wildcards.organelle,
                )
            )
        output:
            (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.pycirclize.pdf"
            )
        log:
            out = "logs/plot_organelle_pycirclize_{organelle}_{tool}_{assembly_name}.out",
            err = "logs/plot_organelle_pycirclize_{organelle}_{tool}_{assembly_name}.err"
        conda:
            "../envs/pycirclize.yml"
        wildcard_constraints:
            organelle = "mitochondrion|chloroplast",
            tool = "pmga|pga_v2"
        params:
            assembly_name = lambda wildcards: wildcards.assembly_name,
            organelle = lambda wildcards: wildcards.organelle,
            label_features = True
        script:
            "../scripts/organelle_pycirclize.py"


if configured_oatk_organelles_with_annotation():
    rule plot_organelle_gbdraw:
        input:
            annotation = (
                lambda wildcards: organelle_annotation_input_for_downstream(
                    wildcards.assembly_name,
                    wildcards.organelle,
                )
            )
        output:
            (
                "results/organelle_annotation/{organelle}/{tool}/{assembly_name}/"
                "{assembly_name}.{organelle}.gbdraw.pdf"
            )
        log:
            out = "logs/plot_organelle_gbdraw_{organelle}_{tool}_{assembly_name}.out",
            err = "logs/plot_organelle_gbdraw_{organelle}_{tool}_{assembly_name}.err"
        conda:
            "../envs/gbdraw.yml"
        wildcard_constraints:
            organelle = "mitochondrion|chloroplast",
            tool = "pmga|pga_v2"
        params:
            output_prefix = lambda wildcards: (
                "results/organelle_annotation/"
                f"{wildcards.organelle}/{wildcards.tool}/{wildcards.assembly_name}/"
                f"{wildcards.assembly_name}.{wildcards.organelle}.gbdraw"
            ),
            species = lambda wildcards: wildcards.assembly_name.replace("_", " "),
            labels = "both"
        shell:
            """
            python3 workflow/scripts/organelle_gbdraw.py \
                --gbk {input.annotation:q} \
                --output-prefix {params.output_prefix:q} \
                --formats pdf \
                --species {params.species:q} \
                --labels {params.labels:q} \
                --overwrite > {log.out:q} 2> {log.err:q}
            """
