import glob
import os
import shlex

RAW_FASTQ_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)


def make_pattern_suffix_pairs(pattern_template, suffix_template, suffixes):
    return [
        (pattern_template.format(suffix=suffix), suffix_template.format(suffix=suffix))
        for suffix in suffixes
    ]


def discover_sample_ids(pattern_suffix_pairs):
    sample_ids = set()
    for pattern, suffix in pattern_suffix_pairs:
        for path in glob.glob(pattern):
            basename = os.path.basename(path)
            if basename.endswith(suffix):
                sample_ids.add(basename[:-len(suffix)])
    return sorted(sample_ids)


def require_sample_ids(sample_ids, description, search_paths):
    if sample_ids:
        return sample_ids
    raise RuntimeError(
        f"No {description} were found. Looked for inputs or intermediates in: {', '.join(search_paths)}"
    )


def resolve_unique_existing_path(candidate_paths, description):
    existing_paths = sorted(path for path in candidate_paths if os.path.exists(path))
    if len(existing_paths) == 1:
        return existing_paths[0]
    if len(existing_paths) > 1:
        raise RuntimeError(
            f"Multiple {description} files were found: {', '.join(existing_paths)}. "
            "Please keep only one matching input file."
        )
    raise RuntimeError(
        f"No {description} file was found. Looked for: {', '.join(candidate_paths)}"
    )


VALID_SELECTED_ASSEMBLIES = (
    "primary",
    "hap1",
    "hap2",
)

HIFIASM_SELECTED_ASSEMBLY_GFA_PATHS = {
    "default": {
        "primary": "results/hifiasm/hifiasm/{assembly_name}.asm.bp.p_ctg.gfa",
        "hap1": "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap1.p_ctg.gfa",
        "hap2": "results/hifiasm/hifiasm/{assembly_name}.asm.bp.hap2.p_ctg.gfa",
    },
    "hic": {
        "primary": "results/hifiasm/hifiasm/{assembly_name}.asm.hic.p_ctg.gfa",
        "hap1": "results/hifiasm/hifiasm/{assembly_name}.asm.hic.hap1.p_ctg.gfa",
        "hap2": "results/hifiasm/hifiasm/{assembly_name}.asm.hic.hap2.p_ctg.gfa",
    },
}


def validate_choice_list(config_key, values, allowed_values):
    if not isinstance(values, list):
        raise ValueError(
            f"'{config_key}' in config.yml must be a YAML list containing one or more values "
            f"from: {', '.join(allowed_values)}"
        )
    if not values:
        raise ValueError(f"'{config_key}' in config.yml must contain at least one value.")

    normalized_values = []
    seen = set()
    for value in values:
        if not isinstance(value, str):
            raise ValueError(
                f"'{config_key}' in config.yml must contain only strings, but found: {value!r}"
            )
        if value not in allowed_values:
            raise ValueError(
                f"Invalid value '{value}' for '{config_key}' in config.yml. Must be one of: "
                f"{', '.join(allowed_values)}"
            )
        if value in seen:
            raise ValueError(
                f"Duplicate value '{value}' found in '{config_key}' in config.yml."
            )
        normalized_values.append(value)
        seen.add(value)
    return normalized_values


def normalize_optional_path_list(config_key, values):
    if values is None:
        return None
    if isinstance(values, str):
        values = [values]
    if not isinstance(values, list):
        raise ValueError(
            f"'{config_key}' in config.yml must be null, a string path, or a YAML list "
            "containing one or more string paths."
        )
    if not values:
        raise ValueError(f"'{config_key}' in config.yml must not be an empty list.")

    normalized_values = []
    for value in values:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                f"'{config_key}' in config.yml must contain only non-empty string paths, "
                f"but found: {value!r}"
            )
        normalized_values.append(value)
    return normalized_values


selected_assemblies = validate_choice_list(
    "selected_assemblies",
    config.get("selected_assemblies", ["primary"]),
    VALID_SELECTED_ASSEMBLIES,
)
submission_assemblies = validate_choice_list(
    "submission_assemblies",
    config.get("submission_assemblies", ["primary"]),
    VALID_SELECTED_ASSEMBLIES,
)
if not set(submission_assemblies).issubset(selected_assemblies):
    raise ValueError(
        "'submission_assemblies' in config.yml must be a subset of 'selected_assemblies'."
    )

selected_assembly_pattern = "|".join(selected_assemblies)
submission_assembly_pattern = "|".join(submission_assemblies)

hic_reads_r1 = normalize_optional_path_list(
    "hic_reads_r1",
    config.get("hic_reads_r1", None),
)
hic_reads_r2 = normalize_optional_path_list(
    "hic_reads_r2",
    config.get("hic_reads_r2", None),
)
if (hic_reads_r1 is None) != (hic_reads_r2 is None):
    raise ValueError(
        "'hic_reads_r1' and 'hic_reads_r2' in config.yml must either both be set or both be null."
    )
if hic_reads_r1 is not None and len(hic_reads_r1) != len(hic_reads_r2):
    raise ValueError(
        "'hic_reads_r1' and 'hic_reads_r2' in config.yml must contain the same number of files."
    )

hic_reads_enabled = hic_reads_r1 is not None
hifiasm_selected_assembly_gfa_paths = HIFIASM_SELECTED_ASSEMBLY_GFA_PATHS[
    "hic" if hic_reads_enabled else "default"
]
downstream_assembly_name = "yahs" if hic_reads_enabled else "fcs"


def expand_selected_assembly_paths(pattern, assembly_name, assemblies=None, **extra_wildcards):
    return expand(
        pattern,
        assembly_name=assembly_name,
        selected_assembly=assemblies or selected_assemblies,
        **extra_wildcards,
    )


def hifiasm_selected_assembly_gfa_path(assembly_name, selected_assembly):
    return hifiasm_selected_assembly_gfa_paths[selected_assembly].format(
        assembly_name=assembly_name
    )


def hifiasm_output_path():
    return {
        f"{selected_assembly}_gfa": path
        for selected_assembly, path in hifiasm_selected_assembly_gfa_paths.items()
    }


def shell_quote_paths(paths):
    return " ".join(shlex.quote(str(path)) for path in paths)


def hic_reads_input(pair):
    if not hic_reads_enabled:
        raise RuntimeError(
            "Hi-C reads were requested, but 'hic_reads_r1' and 'hic_reads_r2' are not set in config.yml."
        )
    if pair == 1:
        return hic_reads_r1
    if pair == 2:
        return hic_reads_r2
    raise ValueError(f"Invalid Hi-C read pair: {pair}. Must be 1 or 2.")


def merged_hic_reads_path(assembly_name, pair):
    return f"results/hic_reads/merged/{assembly_name}_hic_R{pair}.fastq.gz"


def downstream_assembly_path(assembly_name, selected_assembly):
    return (
        f"results/{downstream_assembly_name}/assembly/{selected_assembly}/{assembly_name}.fa"
    )


def seqkit_stats_organelle_path():
    mito_txt = "results/oatk/seqkit/{assembly_name}_mito_seqkit_stats.txt"
    mito_tsv = "results/oatk/seqkit/{assembly_name}_mito_seqkit_stats.tsv"
    pltd_txt = "results/oatk/seqkit/{assembly_name}_pltd_seqkit_stats.txt"
    pltd_tsv = "results/oatk/seqkit/{assembly_name}_pltd_seqkit_stats.tsv"
    if config["oatk_organelle"] == "mito":
        return {"mito_txt": mito_txt, "mito_tsv": mito_tsv}
    if config["oatk_organelle"] == "pltd":
        return {"pltd_txt": pltd_txt, "pltd_tsv": pltd_tsv}
    if config["oatk_organelle"] == "mito_and_pltd":
        return {
            "mito_txt": mito_txt,
            "mito_tsv": mito_tsv,
            "pltd_txt": pltd_txt,
            "pltd_tsv": pltd_tsv,
        }
    raise ValueError(
        "Invalid value for 'oatk_organelle' in config.yml. Must be one of 'mito', "
        "'pltd', or 'mito_and_pltd'."
    )


def seqkit_stats_organelle_real_path(assembly_name):
    return [
        value.format(assembly_name=assembly_name)
        for value in seqkit_stats_organelle_path().values()
    ]


def assembly_all_inputs(assembly_name):
    inputs = [
        f"results/hifi_reads/smudgeplot/{assembly_name}_masked_errors_smu.txt",
        f"results/hifi_reads/genomescope2/{assembly_name}_summary.txt",
        *seqkit_stats_organelle_real_path(assembly_name),
    ]
    for pattern in (
        "results/hifiasm/seqkit/{selected_assembly}/{assembly_name}_seqkit_stats.tsv",
        "results/hifiasm/length/{selected_assembly}/{assembly_name}_length.pdf",
        "results/hifiasm/gc_content/{selected_assembly}/{assembly_name}_gc_content.pdf",
        "results/hifiasm/busco_genome/{selected_assembly}/BUSCO_{assembly_name}.fa",
        "results/hifiasm/merqury/{selected_assembly}/{assembly_name}.merqury.qv",
        "results/hifiasm/depth/{selected_assembly}/{assembly_name}_contig_depth.pdf",
        "results/hifiasm/tidk/{selected_assembly}/{assembly_name}_tidk_find.svg",
        "results/hifiasm/tidk/{selected_assembly}/{assembly_name}_tidk_explore.tsv",
        "results/hifiasm/tidk/{selected_assembly}/{assembly_name}_tidk_search.svg",
    ):
        inputs.extend(expand_selected_assembly_paths(pattern, assembly_name))
    return inputs


def remove_organelle_all_inputs(assembly_name):
    inputs = assembly_all_inputs(assembly_name)
    for pattern in (
        "results/organelle_removal/seqkit/{selected_assembly}/{assembly_name}_seqkit_stats.tsv",
        "results/organelle_removal/length/{selected_assembly}/{assembly_name}_length.pdf",
        "results/organelle_removal/gc_content/{selected_assembly}/{assembly_name}_gc_content.pdf",
        "results/organelle_removal/busco_genome/{selected_assembly}/BUSCO_{assembly_name}.fa",
        "results/organelle_removal/merqury/{selected_assembly}/{assembly_name}.merqury.qv",
        "results/organelle_removal/depth/{selected_assembly}/{assembly_name}_contig_depth.pdf",
        "results/organelle_removal/tidk/{selected_assembly}/{assembly_name}_tidk_find.svg",
        "results/organelle_removal/tidk/{selected_assembly}/{assembly_name}_tidk_explore.tsv",
        "results/organelle_removal/tidk/{selected_assembly}/{assembly_name}_tidk_search.svg",
        "results/hifiasm/map_to_organelle/organelle_contig_depth/{selected_assembly}/{assembly_name}_contig_depth.pdf",
    ):
        inputs.extend(expand_selected_assembly_paths(pattern, assembly_name))
    return inputs


def remove_contamination_all_inputs(assembly_name):
    inputs = remove_organelle_all_inputs(assembly_name)
    for pattern in (
        "results/fcs/seqkit/{selected_assembly}/{assembly_name}_seqkit_stats.tsv",
        "results/fcs/length/{selected_assembly}/{assembly_name}_length.pdf",
        "results/fcs/gc_content/{selected_assembly}/{assembly_name}_gc_content.pdf",
        "results/fcs/busco_genome/{selected_assembly}/BUSCO_{assembly_name}.fa",
        "results/fcs/merqury/{selected_assembly}/{assembly_name}.merqury.qv",
        "results/fcs/depth/{selected_assembly}/{assembly_name}_contig_depth.pdf",
    ):
        inputs.extend(expand_selected_assembly_paths(pattern, assembly_name))
    inputs.append(
        f"results/downloads/tidk/.{assembly_name}_.local_share_tidk_successfully_removed_or_restored.txt"
    )
    return inputs


def scaffold_all_inputs(assembly_name):
    inputs = remove_contamination_all_inputs(assembly_name)
    if not hic_reads_enabled:
        return inputs

    inputs.extend(
        [
            merged_hic_reads_path(assembly_name, 1),
            merged_hic_reads_path(assembly_name, 2),
        ]
    )
    for pattern in (
        "results/yahs/assembly/{selected_assembly}/{assembly_name}.fa",
        "results/yahs/agp/{selected_assembly}/{assembly_name}.agp",
        "results/yahs/alignment/{selected_assembly}/{assembly_name}_hic_to_assembly.bam",
        "results/juicebox/{selected_assembly}/{assembly_name}.hic",
        "results/juicebox/{selected_assembly}/{assembly_name}.assembly",
        "results/juicebox/{selected_assembly}/{assembly_name}.liftover.agp",
    ):
        inputs.extend(expand_selected_assembly_paths(pattern, assembly_name))
    return inputs


def softmask_all_inputs(assembly_name):
    return scaffold_all_inputs(assembly_name) + expand_selected_assembly_paths(
        "results/repeatmasker/{selected_assembly}/{assembly_name}.fa.masked",
        assembly_name,
    )


def gene_prediction_all_inputs(assembly_name, assembly_version):
    inputs = softmask_all_inputs(assembly_name)
    for pattern in (
        "results/isoforms/busco_proteins/{selected_assembly}/BUSCO_{assembly_name}_aa.fa",
        "results/longest_cds/busco_proteins/{selected_assembly}/BUSCO_{assembly_name}_aa.fa",
        "results/isoforms/seqkit/{selected_assembly}/{assembly_name}_seqkit_stats.tsv",
        "results/longest_cds/seqkit/{selected_assembly}/{assembly_name}_seqkit_stats.tsv",
        "results/longest_cds/omark/{selected_assembly}/{assembly_name}_omark",
        "results/longest_cds/{selected_assembly}/{assembly_name}_transcript.fa",
    ):
        inputs.extend(expand_selected_assembly_paths(pattern, assembly_name))
    for pattern in (
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_genome.fa.gz",
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_isoforms.cds.fa.gz",
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_isoforms.gff3.gz",
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_representative.cds.fa.gz",
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_representative.gff3.gz",
        "results/submission/{selected_assembly}/{assembly_name}_{assembly_version}_README.md",
    ):
        inputs.extend(
            expand_selected_assembly_paths(
                pattern,
                assembly_name,
                assemblies=submission_assemblies,
                assembly_version=assembly_version,
            )
        )
    return inputs


def circos_plot_all_inputs(assembly_name, assembly_version):
    return gene_prediction_all_inputs(
        assembly_name, assembly_version
    ) + expand_selected_assembly_paths(
        "results/circos_plot/{selected_assembly}/{assembly_name}_{plot_kind}.pdf",
        assembly_name,
        plot_kind=["circos_plot", "linear_plot"],
    )


def hifi_fastplong_reads(_wildcards):
    sample_ids = require_sample_ids(
        hifi_sample_ids,
        "HiFi read samples",
        [
            "raw_data/*.hifi_reads.bam",
            "results/hifi_reads/fastplong/*_hifi_reads_curated.fastq.gz",
        ],
    )
    return expand(
        "results/hifi_reads/fastplong/{hifi_sample_id}_hifi_reads_curated.fastq.gz",
        hifi_sample_id=sample_ids,
    )


def ont_reads_input():
    ont_reads = config.get("ont_reads", None)
    if ont_reads:
        return ont_reads
    raise RuntimeError(
        "ONT reads were requested, but 'ont_reads' is not set in config.yml and no ONT preprocessing target was provided."
    )


def hifiasm_input(wildcards):
    inputs = {
        "hifi": "results/hifi_reads/merged/{assembly_name}_hifi_reads_curated.fastq.gz".format(
            assembly_name=wildcards.assembly_name
        )
    }
    if config.get("ont_reads", None):
        inputs["ont"] = (
            "results/ont_reads/fastplong/{assembly_name}_ont_reads_curated.fastq.gz".format(
                assembly_name=wildcards.assembly_name
            )
        )
    if hic_reads_enabled:
        inputs["hic_r1"] = merged_hic_reads_path(wildcards.assembly_name, 1)
        inputs["hic_r2"] = merged_hic_reads_path(wildcards.assembly_name, 2)
    return inputs


def hifiasm_hic_option(input):
    if "hic_r1" not in input.keys():
        return ""
    return f"--h1 {input.hic_r1} --h2 {input.hic_r2}"


def oatkdb_path():
    mito_fam = f"results/downloads/oatkdb/{config['oatk_lineage']}_mito.fam"
    pltd_fam = f"results/downloads/oatkdb/{config['oatk_lineage']}_pltd.fam"
    if config["oatk_organelle"] == "mito":
        return {"mito_fam": mito_fam}
    if config["oatk_organelle"] == "pltd":
        return {"pltd_fam": pltd_fam}
    if config["oatk_organelle"] == "mito_and_pltd":
        return {"mito_fam": mito_fam, "pltd_fam": pltd_fam}
    raise ValueError(
        "Invalid value for 'oatk_organelle' in config.yml. Must be one of 'mito', "
        "'pltd', or 'mito_and_pltd'."
    )


def oatk_output_path():
    utg_final_gfa = "results/oatk/oatk/{assembly_name}.utg.final.gfa"
    annot_mito_txt = "results/oatk/oatk/{assembly_name}.annot_mito.txt"
    annot_pltd_txt = "results/oatk/oatk/{assembly_name}.annot_pltd.txt"
    mito_gfa = "results/oatk/oatk/{assembly_name}.mito.gfa"
    mito_bed = "results/oatk/oatk/{assembly_name}.mito.bed"
    mito_ctg_fasta = "results/oatk/oatk/{assembly_name}.mito.ctg.fasta"
    mito_ctg_bed = "results/oatk/oatk/{assembly_name}.mito.ctg.bed"
    pltd_gfa = "results/oatk/oatk/{assembly_name}.pltd.gfa"
    pltd_bed = "results/oatk/oatk/{assembly_name}.pltd.bed"
    pltd_ctg_fasta = "results/oatk/oatk/{assembly_name}.pltd.ctg.fasta"
    pltd_ctg_bed = "results/oatk/oatk/{assembly_name}.pltd.ctg.bed"
    if config["oatk_organelle"] == "mito":
        return {
            "utg_final_gfa": utg_final_gfa,
            "annot_mito_txt": annot_mito_txt,
            "mito_gfa": mito_gfa,
            "mito_bed": mito_bed,
            "mito_ctg_fasta": mito_ctg_fasta,
            "mito_ctg_bed": mito_ctg_bed,
        }
    if config["oatk_organelle"] == "pltd":
        return {
            "utg_final_gfa": utg_final_gfa,
            "annot_pltd_txt": annot_pltd_txt,
            "pltd_gfa": pltd_gfa,
            "pltd_bed": pltd_bed,
            "pltd_ctg_fasta": pltd_ctg_fasta,
            "pltd_ctg_bed": pltd_ctg_bed,
        }
    if config["oatk_organelle"] == "mito_and_pltd":
        return {
            "utg_final_gfa": utg_final_gfa,
            "annot_mito_txt": annot_mito_txt,
            "annot_pltd_txt": annot_pltd_txt,
            "mito_gfa": mito_gfa,
            "mito_bed": mito_bed,
            "mito_ctg_fasta": mito_ctg_fasta,
            "mito_ctg_bed": mito_ctg_bed,
            "pltd_gfa": pltd_gfa,
            "pltd_bed": pltd_bed,
            "pltd_ctg_fasta": pltd_ctg_fasta,
            "pltd_ctg_bed": pltd_ctg_bed,
        }
    raise ValueError(
        "Invalid value for 'oatk_organelle' in config.yml. Must be one of 'mito', "
        "'pltd', or 'mito_and_pltd'."
    )


def concatemer_path():
    mito = "results/oatk/concatemer/{assembly_name}.concatemer.mito.fa"
    pltd = "results/oatk/concatemer/{assembly_name}.concatemer.pltd.fa"
    all_organelle = "results/oatk/concatemer/{assembly_name}.concatemer.all.fa"
    if config["oatk_organelle"] == "mito":
        return {"mito": mito, "all_organelle": all_organelle}
    if config["oatk_organelle"] == "pltd":
        return {"pltd": pltd, "all_organelle": all_organelle}
    if config["oatk_organelle"] == "mito_and_pltd":
        return {"mito": mito, "pltd": pltd, "all_organelle": all_organelle}
    raise ValueError(
        "Invalid value for 'oatk_organelle' in config.yml. Must be one of 'mito', "
        "'pltd', or 'mito_and_pltd'."
    )


def rnaseq_raw_input_path(rnaseq_sample_id, pair):
    return resolve_unique_existing_path(
        [
            f"raw_data/{rnaseq_sample_id}_{pair}{suffix}"
            for suffix in RAW_FASTQ_SUFFIXES
        ],
        f"RNA-Seq read pair {pair} for sample '{rnaseq_sample_id}'",
    )


def braker3_rnaseq_inputs(_wildcards):
    sample_ids = require_sample_ids(
        rnaseq_sample_ids,
        "RNA-Seq samples",
        [
            *(pattern for pattern, _ in rnaseq_raw_pattern_suffix_pairs),
            rnaseq_fastp_search_path,
        ],
    )
    return expand(
        "results/rnaseq_reads/fastp/{rnaseq_sample_id}_{pair}.fastq",
        rnaseq_sample_id=sample_ids,
        pair=[1, 2],
    )


def circos_plot_input_path(wildcards):
    selected_assembly = wildcards.selected_assembly
    assembly_name = wildcards.assembly_name
    input_dict = {
        "contig": (
            f"results/circos_plot/{selected_assembly}/{assembly_name}_long_contig_length.tsv"
        )
    }
    if "gene" in circos_ids:
        input_dict["gene"] = (
            f"results/circos_plot/gene/{selected_assembly}/{assembly_name}_windows_gene_coverage.bed"
        )
    for repeat_class in repeat_classes:
        input_dict[repeat_class] = (
            f"results/circos_plot/repeatmasker_repeat/{selected_assembly}/{repeat_class}/"
            f"{assembly_name}_windows_repeat_coverage.bed"
        )
    if "tidk" in circos_ids:
        input_dict["tidk"] = (
            f"results/circos_plot/tidk_repeat/{selected_assembly}/{assembly_name}_windows_tidk.tsv"
        )
    return input_dict
