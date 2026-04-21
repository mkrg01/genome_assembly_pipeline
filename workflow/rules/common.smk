import glob
import os

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
    return [
        f"results/hifi_reads/smudgeplot/{assembly_name}_masked_errors_smu.txt",
        f"results/hifi_reads/genomescope2/{assembly_name}_summary.txt",
        *seqkit_stats_organelle_real_path(assembly_name),
        f"results/hifiasm/seqkit/{assembly_name}_seqkit_stats.tsv",
        f"results/hifiasm/length/{assembly_name}_length.pdf",
        f"results/hifiasm/gc_content/{assembly_name}_gc_content.pdf",
        f"results/hifiasm/busco_genome/BUSCO_{assembly_name}.asm.bp.p_ctg.fa",
        f"results/hifiasm/merqury/{assembly_name}.merqury.qv",
        f"results/hifiasm/depth/{assembly_name}_contig_depth.pdf",
        f"results/hifiasm/tidk/{assembly_name}_tidk_find.svg",
        f"results/hifiasm/tidk/{assembly_name}_tidk_explore.tsv",
        f"results/hifiasm/tidk/{assembly_name}_tidk_search.svg",
    ]


def remove_organelle_all_inputs(assembly_name):
    return assembly_all_inputs(assembly_name) + [
        f"results/organelle_removal/seqkit/{assembly_name}_seqkit_stats.tsv",
        f"results/organelle_removal/length/{assembly_name}_length.pdf",
        f"results/organelle_removal/gc_content/{assembly_name}_gc_content.pdf",
        f"results/organelle_removal/busco_genome/BUSCO_{assembly_name}.asm.bp.p_ctg.fa",
        f"results/organelle_removal/merqury/{assembly_name}.merqury.qv",
        f"results/organelle_removal/depth/{assembly_name}_contig_depth.pdf",
        f"results/organelle_removal/tidk/{assembly_name}_tidk_find.svg",
        f"results/organelle_removal/tidk/{assembly_name}_tidk_explore.tsv",
        f"results/organelle_removal/tidk/{assembly_name}_tidk_search.svg",
        f"results/hifiasm/map_to_organelle/organelle_contig_depth/{assembly_name}_contig_depth.pdf",
    ]


def remove_contamination_all_inputs(assembly_name):
    return remove_organelle_all_inputs(assembly_name) + [
        f"results/fcs/seqkit/{assembly_name}_seqkit_stats.tsv",
        f"results/fcs/length/{assembly_name}_length.pdf",
        f"results/fcs/gc_content/{assembly_name}_gc_content.pdf",
        f"results/fcs/busco_genome/BUSCO_{assembly_name}.asm.bp.p_ctg.fa",
        f"results/fcs/merqury/{assembly_name}.merqury.qv",
        f"results/fcs/depth/{assembly_name}_contig_depth.pdf",
        f"results/downloads/tidk/.{assembly_name}_.local_share_tidk_successfully_removed_or_restored.txt",
    ]


def softmask_all_inputs(assembly_name):
    return remove_contamination_all_inputs(assembly_name) + [
        f"results/repeatmasker/{assembly_name}.asm.bp.p_ctg.fa.masked",
    ]


def gene_prediction_all_inputs(assembly_name, assembly_version):
    return softmask_all_inputs(assembly_name) + [
        f"results/isoforms/busco_proteins/BUSCO_{assembly_name}_aa.fa",
        f"results/longest_cds/busco_proteins/BUSCO_{assembly_name}_aa.fa",
        f"results/isoforms/seqkit/{assembly_name}_seqkit_stats.tsv",
        f"results/longest_cds/seqkit/{assembly_name}_seqkit_stats.tsv",
        f"results/longest_cds/omark/{assembly_name}_omark",
        f"results/longest_cds/{assembly_name}_transcript.fa",
        f"results/submission/{assembly_name}_{assembly_version}_genome.fa.gz",
        f"results/submission/{assembly_name}_{assembly_version}_isoforms.cds.fa.gz",
        f"results/submission/{assembly_name}_{assembly_version}_isoforms.gff3.gz",
        f"results/submission/{assembly_name}_{assembly_version}_representative.cds.fa.gz",
        f"results/submission/{assembly_name}_{assembly_version}_representative.gff3.gz",
        f"results/submission/{assembly_name}_{assembly_version}_README.md",
    ]


def circos_plot_all_inputs(assembly_name, assembly_version):
    return gene_prediction_all_inputs(assembly_name, assembly_version) + [
        f"results/circos_plot/{assembly_name}_circos_plot.pdf",
        f"results/circos_plot/{assembly_name}_linear_plot.pdf",
    ]


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
    return inputs


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


def circos_plot_input_path():
    input_dict = {
        "contig": "results/circos_plot/{assembly_name}_long_contig_length.tsv"
    }
    if "gene" in circos_ids:
        input_dict["gene"] = (
            "results/circos_plot/gene/{assembly_name}_windows_gene_coverage.bed"
        )
    for repeat_class in repeat_classes:
        input_dict[repeat_class] = (
            f"results/circos_plot/repeatmasker_repeat/{repeat_class}/"
            "{assembly_name}_windows_repeat_coverage.bed"
        )
    if "tidk" in circos_ids:
        input_dict["tidk"] = (
            "results/circos_plot/tidk_repeat/{assembly_name}_windows_tidk.tsv"
        )
    return input_dict
