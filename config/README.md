# Configuration Guide

This document explains:

1. **Required input files** to be placed in the `raw_data` directory.  
2. **Configuration parameters** to set in `config/config.yml`.

---

## 1. Input Files (`raw_data/`)

The full assembly workflow requires **PacBio HiFi reads**. The external-assembly annotation workflow (`workflow/Snakefile.annotation`) instead requires an assembly FASTA via `external_assembly`. **Paired-end RNA-seq reads** are required when running targets that include gene prediction (`gene_prediction_all`, `circos_plot_all`, or `all`).

Place your raw sequencing files in the `raw_data` directory with the following naming conventions:

| File Type                  | Naming Pattern               | Example                                    |
| -------------------------- | ---------------------------- | ------------------------------------------ |
| PacBio HiFi reads          | `*.hifi_reads.bam` | `SAMPLE1.hifi_reads.bam`                    |
| Index for HiFi reads  | `*.hifi_reads.bam.pbi` | `SAMPLE1.hifi_reads.bam.pbi`                    |
| External assembly FASTA (annotation workflow only) | Any path set in `external_assembly` | `raw_data/GCA_054852875.1_ASM5485287v1_genomic.fna` |
| Ultra-long ONT reads (optional) | Any FASTQ path in config (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/ont_reads.fastq.gz`              |
| Hi-C reads R1 (optional) | One or more FASTQ paths in `hic_reads_r1` (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/hic_R1.fastq.gz` |
| Hi-C reads R2 (optional) | One or more FASTQ paths in `hic_reads_r2` (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) | `raw_data/hic_R2.fastq.gz` |
| Paired-end RNA-seq (R1)    | `*_1.fastq.gz`, `*_1.fq.gz`, `*_1.fastq`, or `*_1.fq`        | `RNASEQ1_1.fastq.gz`                       |
| Paired-end RNA-seq (R2)    | `*_2.fastq.gz`, `*_2.fq.gz`, `*_2.fastq`, or `*_2.fq`        | `RNASEQ1_2.fastq.gz`                       |

**Notes:**

- The pipeline will automatically detect and process multiple PacBio HiFi and RNA-seq samples, if present.
- PacBio HiFi reads are not required when running the external-assembly annotation workflow with `workflow/Snakefile.annotation`.
- RNA-seq raw FASTQ inputs can be gzipped or plain text. The pipeline normalizes them during preprocessing, so downstream rules keep using the same internal filenames.
- Keep only one raw RNA-seq file per sample and read pair. For example, do not place both `SAMPLE_1.fastq.gz` and `SAMPLE_1.fq` in `raw_data/`.
- Ultra-long ONT reads are optional and can improve assembly quality when integrated with HiFi reads. See the [hifiasm documentation](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#ultra-long-ont-integration) for details.
- Hi-C reads are optional and should be provided as matching R1/R2 paths in `hic_reads_r1` and `hic_reads_r2`. Multiple lanes can be listed in matching order. If Hi-C reads are configured, the pipeline uses them for hifiasm phasing, runs YaHS scaffolding on each assembly listed in `selected_assemblies`, and prepares Juicebox-ready contact maps. See the [hifiasm documentation](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#hi-c-integration) and the [YaHS documentation](https://github.com/c-zhou/yahs) for details.
- LongStitch, when enabled, uses the curated HiFi reads generated internally by this workflow.

---

## 2. Configuration File (`config/config.yml`)

Edit `config/config.yml` to match your dataset and analysis requirements.  
The sections below follow the order of `config/config.yml`.

### Runtime

| Parameter | Description | Example |
| --- | --- | --- |
| `pipeline_version` | Container image tag used by the entrypoint Snakefiles. The workflow builds the full image URI as `docker://aurelia01/genome_assembly_pipeline:<pipeline_version>`. | `"v0.6.15"` |
| `organism_name` | Filesystem-safe organism label used for output file prefixes and GenBank source/metadata. Use underscores instead of spaces; legacy `assembly_name` is still accepted as a fallback. | `"Dioncophyllum_thollonii"` |
| `genome_version` | Filesystem-safe version label for this genome release. Used for release directories and output file names; legacy `assembly_version` is still accepted as a fallback. | `"v1.0"` |
| `external_assembly` | Required only with `workflow/Snakefile.annotation`. Set this to a single FASTA path string. Gzipped FASTA inputs are decompressed when staged into `results/external/assembly/`. | `"raw_data/GCA_054852875.1_ASM5485287v1_genomic.fna"` |

### Hifiasm Genome Assembly

| Parameter | Description | Example |
| --- | --- | --- |
| `selected_assemblies`   | One or more assemblies to process after hifiasm. Allowed values: `{"primary", "hap1", "hap2"}`. Downstream outputs are written under assembly-specific subdirectories such as `results/fcs/assembly/<selected_assembly>/...` or `results/release/<genome_version>/nuclear/<selected_assembly>/...`. Release outputs are produced for every selected assembly. When Hi-C integration is enabled, these correspond to hifiasm `hic.p_ctg`, `hic.hap1.p_ctg`, and `hic.hap2.p_ctg` outputs, which are then scaffolded with YaHS, exported as Juicebox-ready contact-map files under `results/juicebox/<selected_assembly>/...`, and written under `results/yahs/assembly/<selected_assembly>/...` before RepeatMasker and later steps. Organelle release assets are staged separately under `results/release/<genome_version>/organelle/...` according to `oatk_organelle` and `organelle_annotation`. | `["primary", "hap1", "hap2"]` |
| `hifiasm_dual_scaf`    | Optional: Enable hifiasm self-scaffolding by adding `--dual-scaf`. `{true, false}` | `false` |
| `ont_reads`             | Optional: Path to ultra-long ONT reads in FASTQ format (`.fastq.gz`, `.fq.gz`, `.fastq`, or `.fq`). Set to `null` to disable ONT integration. [See hifiasm docs](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#ultra-long-ont-integration) | `null` |
| `hic_reads_r1`          | Optional: Hi-C read 1 input for hifiasm phasing, YaHS scaffolding, and Juicebox-ready contact maps. Set to `null` to disable Hi-C integration, or provide one or more FASTQ paths as a YAML list. Must be paired with `hic_reads_r2`. [See hifiasm docs](https://github.com/chhylp123/hifiasm?tab=readme-ov-file#hi-c-integration) | `["raw_data/hic_R1.fastq.gz"]` |
| `hic_reads_r2`          | Optional: Hi-C read 2 input for hifiasm phasing, YaHS scaffolding, and Juicebox-ready contact maps. Must match `hic_reads_r1` in length and order. | `["raw_data/hic_R2.fastq.gz"]` |

### Oatk Organelle Assembly

| Parameter | Description | Example |
| --- | --- | --- |
| `oatk_lineage`          | Lineage of the Oatk HMM profile database. [Lineage list](https://github.com/c-zhou/OatkDB/blob/main/v20230921/TAXID) | `"magnoliopsida"` |
| `oatk_organelle`        | Organelle to assemble. `{"mitochondrion", "chloroplast", "mitochondrion_and_chloroplast"}`. Legacy aliases `mito`, `pltd`, and `mito_and_pltd` are still accepted. | `"mitochondrion_and_chloroplast"` |
| `oatk_minimum_kmer_coverage`| Minimum kmer coverage used for Oatk. [Instructions](https://github.com/c-zhou/oatk)  | `"250"` |

### Contamination Screening

| Parameter | Description | Example |
| --- | --- | --- |
| `taxid`          | NCBI Taxonomy ID for the target organism. Used by FCS-GX screening and organelle GenBank source `db_xref`. Legacy `fcs_gx_taxid` is still accepted as a fallback. [NCBI Taxonomy Tree](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/) | `"122299"` for *Dioncophyllum thollonii* |

### LongStitch Misassembly Correction

| Parameter | Description | Example |
| --- | --- | --- |
| `longstitch_enabled` | Optional: Run LongStitch on the FCS-cleaned assembly before YaHS and downstream analysis. LongStitch uses the curated HiFi reads generated by this workflow and always runs with the HiFi long-read mapping preset. When enabled, the haploid genome size passed to LongStitch as `G` is derived from the FCS-cleaned assembly `sum_len` in the workflow-generated `seqkit stats` table. `{true, false}` | `false` |

### YaHS Scaffolding

| Parameter | Description | Example |
| --- | --- | --- |
| `yahs_restriction_enzymes` | Optional: Restriction enzyme motif(s) passed to YaHS with `-e` during Hi-C scaffolding after FCS cleanup and optional LongStitch correction. Set to `null` to use YaHS defaults. This parameter is used only when both `hic_reads_r1` and `hic_reads_r2` are set. [See YaHS docs](https://github.com/c-zhou/yahs) | `"GATC"` |

### Quality Assessment

| Parameter | Description | Example |
| --- | --- | --- |
| `busco_lineage_dataset` | BUSCO lineage dataset for genome completeness assessment. [Lineage list](https://busco-data.ezlab.org/v5/data/lineages/) | `"embryophyta_odb12"`                      |
| `tidk_clade`            | Clade for [tidk find](https://github.com/tolkit/telomeric-identifier). [Lineage list](https://github.com/tolkit/telomeric-identifier?tab=readme-ov-file#find) | `"Caryophyllales"`|
| `tidk_telomeric_repeat_unit` | A telomeric repeat unit for [tidk search](https://github.com/tolkit/telomeric-identifier). [A Telomeric Repeat Database](https://github.com/tolkit/a-telomeric-repeat-database) | `"AAACCCT"`|

### Softmasking

| Parameter | Description | Example |
| --- | --- | --- |
| `dfam_version`          | Version of the Dfam database for RepeatMasker. [Dfam releases](https://www.dfam.org) | `"3.9"`                                    |
| `dfam_partitions`       | Dfam partitions. See [README.txt](https://www.dfam.org/releases/current/families/FamDB/README.txt). | `"0,5,6"` (Viridiplantae)                      |
| `dfam_lineage_name`     | Name of the Dfam lineage to use.                             | `"Viridiplantae"`                          |

### Gene Prediction

| Parameter | Description | Example |
| --- | --- | --- |
| `orthodb_version`       | Version of the OrthoDB database (used by Braker3). [ProtHint instructions](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) | `"12"`                                     |
| `orthodb_lineage`       | OrthoDB lineage dataset to use. [Lineage list](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"Viridiplantae"`                          |
| `orthodb_md5sum`        | MD5 checksum of the OrthoDB database. [Checksums](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"34c1f027a1a7b10f225b69fbd5500587"`       |

### Visualization

| Parameter | Description | Example |
| --- | --- | --- |
| `min_long_contig_length` | Minimum contig length to be considered a "long contig" for downstream visualization (e.g., Circos plots). | `1_000_000` |
| `circos_plot_tracks`        | Track configuration for the Circos plot | `{id: "gene", label: "Gene model", color: "#4C72B0", window_size: 500_000}`       |
| `circos_plot_x_major_tick_interval` | Major tick interval (in bp) for the x-axis in Circos plots. Major ticks are labeled. | `10_000_000` |
| `circos_plot_x_minor_tick_interval` | Minor tick interval (in bp) for the x-axis in Circos plots. Minor ticks are unlabeled. | `5_000_000` |

### Organelle Annotation

| Parameter | Description | Example |
| --- | --- | --- |
| `organelle_annotation` | Annotation tools for Oatk-assembled organelle genomes. Mitochondrion: `{"pmga", null}`. Chloroplast: `{"pga_v2", null}`. Omitted keys, per-organelle `null`, or top-level `organelle_annotation: null` skip annotation. Legacy keys `mito` and `pltd` are still accepted. Oatk organelle FASTA record IDs are prefixed before annotation/submission (`mt_` for mitochondrion, `cp_` for chloroplast) so GenBank LOCUS names and submission FASTA headers stay unique and consistent. For PMGA/PGA, circular single-record inputs are first annotated once to choose an internal feature-free origin, then rotated before final annotation so origin-spanning genes are less likely to be missed; the rotated annotation FASTA is used for RNA-editing evidence mapping and submission. PMGA/PGA GenBank outputs are then post-curated with RNA-seq and HiFi evidence to add CDS `/translation` and accepted `/exception="RNA editing"` qualifiers. Normal RNA-editing calls use accept thresholds of RNA depth >= 10, edited reads >= 3, edit fraction >= 0.10, base quality >= 30, mapping quality >= 30, DNA/HiFi depth >= 10, and DNA/HiFi alternate fraction <= 0.10; CDS-essential rescue sites that restore a valid start codon, terminal stop codon, or premature-stop rescue may use a separate one-read rescue threshold with DNA/HiFi alternate fraction <= 0.10. Sites with DNA/HiFi alternate support above the threshold are retained as `likely_genomic_variant` evidence and summarized for genome sequence or variant review. PMGA and PGA v2.0 are downloaded automatically only when selected. | `{mitochondrion: "pmga", chloroplast: "pga_v2"}` |
| `organelle_reference_cds_qc` | Optional reference CDS/protein QC for PMGA/PGA annotations. Set an organelle `reference_dir` to `null` to skip reference QC and reference-inferred RNA-editing rescue for that organelle. When `reference_dir` is set, the pipeline emits pre/post RNA-editing QC tables and conservatively tries closely related reference CDS/protein evidence for residual invalid CDS after RNA-seq curation. `fix_hifi_frameshifts: true` is allowed only for chloroplast/PGA and may correct well-supported HiFi read-supported frameshift indels before rerunning PGA. Chloroplast/PGA requires a reference directory because PGA v2 itself is reference-based. | `{chloroplast: {reference_dir: "plastid_reference", fix_hifi_frameshifts: false}, mitochondrion: {reference_dir: null, fix_hifi_frameshifts: false}}` |

#### Choosing Organelle Reference CDS QC References

When `organelle_annotation.chloroplast` is set to `pga_v2`, place one or a few plastid GenBank files (`.gb`, `.gbk`, or `.gbff`) in `organelle_reference_cds_qc.chloroplast.reference_dir`. PGA v2 uses these annotated references to transfer chloroplast/plastid gene annotations, so the reference choice can affect annotation completeness and naming. The same references are also used for chloroplast CDS/protein QC and reference-inferred RNA-editing rescue.

For PMGA/mitochondrion, `organelle_reference_cds_qc.mitochondrion.reference_dir` is optional. If provided, place one or a few mitochondrial GenBank files (`.gb`, `.gbk`, or `.gbff`) in that directory. Those references are used for mitochondrial CDS/protein QC and reference-inferred RNA-editing rescue. Applied reference-inferred sites are recorded separately from RNA-seq-supported calls in the RNA-editing evidence sidecars and GenBank `/inference` qualifiers.

Useful NCBI Nucleotide searches for chloroplast/plastid references:

```text
"<Species name>"[Organism] AND chloroplast[Title] AND "complete genome"[Title]
<Genus name>[Organism] AND chloroplast[Title] AND "complete genome"[Title]
<Family name>[Organism] AND chloroplast[Title] AND "complete genome"[Title]
```

Useful NCBI Nucleotide searches for mitochondrial references:

```text
"<Species name>"[Organism] AND mitochondrion[Title] AND "complete genome"[Title]
<Genus name>[Organism] AND mitochondrion[Title] AND "complete genome"[Title]
<Family name>[Organism] AND mitochondrion[Title] AND "complete genome"[Title]
```

Prefer complete organelle GenBank records from the same species, genus, or family when available.

NCBI GenBank download examples:

```bash
mkdir -p plastid_reference
mkdir -p mitochondrial_reference

PLASTID_ACCESSION="NC_041245.1"
MITOCHONDRION_ACCESSION="<mitochondrion accession>"

curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${PLASTID_ACCESSION}&rettype=gb&retmode=text" \
  -o "plastid_reference/${PLASTID_ACCESSION}.gbk"

curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${MITOCHONDRION_ACCESSION}&rettype=gb&retmode=text" \
  -o "mitochondrial_reference/${MITOCHONDRION_ACCESSION}.gbk"
```
