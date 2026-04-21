# Configuration Guide

This document explains:

1. **Required input files** to be placed in the `raw_data` directory.  
2. **Configuration parameters** to set in `config/config.yml`.

---

## 1. Input Files (`raw_data/`)

The pipeline requires **PacBio HiFi reads**. **Paired-end RNA-seq reads** are required when running targets that include gene prediction (`gene_prediction_all`, `circos_plot_all`, or `all`).

Place your raw sequencing files in the `raw_data` directory with the following naming conventions:

| File Type                  | Naming Pattern               | Example                                    |
| -------------------------- | ---------------------------- | ------------------------------------------ |
| PacBio HiFi reads          | `*.hifi_reads.bam` | `SAMPLE1.hifi_reads.bam`                    |
| Index for HiFi reads  | `*.hifi_reads.bam.pbi` | `SAMPLE1.hifi_reads.bam.pbi`                    |
| Ultra-long ONT reads (optional) | Any FASTQ path in config (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/ont_reads.fastq.gz`              |
| Hi-C reads R1 (optional) | One or more FASTQ paths in `hic_reads_r1` (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/hic_R1.fastq.gz` |
| Hi-C reads R2 (optional) | One or more FASTQ paths in `hic_reads_r2` (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/hic_R2.fastq.gz` |
| Paired-end RNA-seq (R1)    | `*_1.fastq.gz`, `*_1.fq.gz`, `*_1.fastq`, or `*_1.fq`        | `RNASEQ1_1.fastq.gz`                       |
| Paired-end RNA-seq (R2)    | `*_2.fastq.gz`, `*_2.fq.gz`, `*_2.fastq`, or `*_2.fq`        | `RNASEQ1_2.fastq.gz`                       |

**Notes:**

- The pipeline will automatically detect and process multiple PacBio HiFi and RNA-seq samples, if present.
- RNA-seq raw FASTQ inputs can be gzipped or plain text. The pipeline normalizes them during preprocessing, so downstream rules keep using the same internal filenames.
- Keep only one raw RNA-seq file per sample and read pair. For example, do not place both `SAMPLE_1.fastq.gz` and `SAMPLE_1.fq` in `raw_data/`.
- Ultra-long ONT reads are optional and can improve assembly quality when integrated with HiFi reads. See the [hifiasm documentation](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#ultra-long-ont-integration) for details.
- Hi-C reads are optional and should be provided as matching R1/R2 paths in `hic_reads_r1` and `hic_reads_r2`. Multiple lanes can be listed in matching order. If Hi-C reads are configured, the pipeline uses them for hifiasm phasing, runs YaHS scaffolding on each assembly listed in `selected_assemblies`, and prepares Juicebox-ready contact maps. See the [hifiasm documentation](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#hi-c-integration) and the [YaHS documentation](https://github.com/c-zhou/yahs) for details.

---

## 2. Configuration File (`config/config.yml`)

Edit `config/config.yml` to match your dataset and analysis requirements.  
Below are the available parameters:

| Parameter               | Description                                                  | Example                                    |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------ |
| `assembly_name`         | Name used for output files | `"Dioncophyllum_thollonii"` |
| `assembly_version`         | Version used for output files | `"v1.0"` |
| `selected_assemblies`   | One or more assemblies to process after hifiasm. Allowed values: `{"primary", "hap1", "hap2"}`. Downstream outputs are written under assembly-specific subdirectories such as `results/fcs/assembly/<selected_assembly>/...` or `results/submission/<selected_assembly>/...`. When Hi-C integration is enabled, these correspond to hifiasm `hic.p_ctg`, `hic.hap1.p_ctg`, and `hic.hap2.p_ctg` outputs, which are then scaffolded with YaHS, exported as Juicebox-ready contact-map files under `results/juicebox/<selected_assembly>/...`, and written under `results/yahs/assembly/<selected_assembly>/...` before RepeatMasker and later steps. | `["primary", "hap1", "hap2"]` |
| `submission_assemblies` | One or more assemblies to package for submission. Must be a subset of `selected_assemblies`. | `["primary"]` |
| `hifiasm_dual_scaf`    | Optional: Enable hifiasm self-scaffolding by adding `--dual-scaf`. The hifiasm documentation describes this as improving contiguity for diploid haplotype-resolved assemblies by letting the two haplotypes scaffold each other. Keep `false` if you want the current default behavior or if self-scaffolding is not appropriate for your sample. | `false` |
| `ont_reads`             | Optional: Path to ultra-long ONT reads in FASTQ format (`.fastq.gz`, `.fq.gz`, `.fastq`, or `.fq`). Set to `null` to disable ONT integration. [See hifiasm docs](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#ultra-long-ont-integration) | `null` |
| `hic_reads_r1`          | Optional: Hi-C read 1 input for hifiasm phasing, YaHS scaffolding, and Juicebox-ready contact maps. Set to `null` to disable Hi-C integration, or provide one or more FASTQ paths as a YAML list. Must be paired with `hic_reads_r2`. [See hifiasm docs](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#hi-c-integration) | `["raw_data/hic_R1.fastq.gz"]` |
| `hic_reads_r2`          | Optional: Hi-C read 2 input for hifiasm phasing, YaHS scaffolding, and Juicebox-ready contact maps. Must match `hic_reads_r1` in length and order. | `["raw_data/hic_R2.fastq.gz"]` |
| `yahs_restriction_enzymes` | Optional: Restriction enzyme motif(s) passed to YaHS with `-e` during Hi-C scaffolding. Set to `null` to use YaHS defaults. [See YaHS docs](https://github.com/c-zhou/yahs) | `"GATC"` |
| `oatk_lineage`          | Lineage of the Oatk HMM profile database. [Lineage list](https://github.com/c-zhou/OatkDB/blob/main/v20230921/TAXID) | `"magnoliopsida"` |
| `oatk_organelle`        | Organelle to assemble. `{"mito", "pltd", "mito_and_pltd"}`  | `"mito_and_pltd"` |
| `oatk_minimum_kmer_coverage`| Minimum kmer coverage used for Oatk. [Instructions](https://github.com/c-zhou/oatk)  | `"250"` |
| `fcs_gx_taxid`          | NCBI Taxonomy ID for FCS-GX screening. [NCBI Taxonomy Tree](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/) | `"122299"` for *Dioncophyllum thollonii* |
| `busco_lineage_dataset` | BUSCO lineage dataset for genome completeness assessment. [Lineage list](https://busco-data.ezlab.org/v5/data/lineages/) | `"embryophyta_odb12"`                      |
| `tidk_clade`            | Clade for [tidk find](https://github.com/tolkit/telomeric-identifier). [Lineage list](https://github.com/tolkit/telomeric-identifier?tab=readme-ov-file#find) | `"Caryophyllales"`|
| `tidk_telomeric_repeat_unit` | A telomeric repeat unit for [tidk search](https://github.com/tolkit/telomeric-identifier). [A Telomeric Repeat Database](https://github.com/tolkit/a-telomeric-repeat-database) | `"AAACCCT"`|
| `dfam_version`          | Version of the Dfam database for RepeatMasker. [Dfam releases](https://www.dfam.org) | `"3.9"`                                    |
| `dfam_partitions`       | Dfam partitions. See [README.txt](https://www.dfam.org/releases/current/families/FamDB/README.txt). | `"0,5,6"` (Viridiplantae)                      |
| `dfam_lineage_name`     | Name of the Dfam lineage to use.                             | `"Viridiplantae"`                          |
| `orthodb_version`       | Version of the OrthoDB database (used by Braker3). [ProtHint instructions](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) | `"12"`                                     |
| `orthodb_lineage`       | OrthoDB lineage dataset to use. [Lineage list](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"Viridiplantae"`                          |
| `orthodb_md5sum`        | MD5 checksum of the OrthoDB database. [Checksums](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"34c1f027a1a7b10f225b69fbd5500587"`       |
| `min_long_contig_length` | Minimum contig length to be considered a "long contig" for downstream visualization (e.g., Circos plots). | `1_000_000` |
| `circos_plot_tracks`        | Track configuration for the Circos plot | `{id: "gene", label: "Gene model", color: "#4C72B0", window_size: 500_000}`       |
| `circos_plot_x_major_tick_interval` | Major tick interval (in bp) for the x-axis in Circos plots. Major ticks are labeled. | `10_000_000` |
| `circos_plot_x_minor_tick_interval` | Minor tick interval (in bp) for the x-axis in Circos plots. Minor ticks are unlabeled. | `5_000_000` |
