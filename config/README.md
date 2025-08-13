# Configuration Guide

This document explains:

1. **Required input files** to be placed in the `raw_data` directory.  
2. **Configuration parameters** to set in `config/config.yaml`.

---

## 1. Input Files (`raw_data/`)

The pipeline requires **both PacBio HiFi reads** and **paired-end RNA-seq reads**.

Place your raw sequencing files in the `raw_data` directory with the following naming conventions:

| File Type                  | Naming Pattern               | Example                                    |
| -------------------------- | ---------------------------- | ------------------------------------------ |
| Paired-end RNA-seq (R1)    | `<RNA_ID>_1.fastq.gz`        | `RNASEQ1_1.fastq.gz`                       |
| Paired-end RNA-seq (R2)    | `<RNA_ID>_2.fastq.gz`        | `RNASEQ1_2.fastq.gz`                       |
| Additional RNA-seq samples | Continue numbering           | `RNASEQ2_1.fastq.gz`, `RNASEQ2_2.fastq.gz` |
| PacBio HiFi reads          | `<SAMPLE_ID>.hifi_reads.bam` | `SAMPLE.hifi_reads.bam`                    |

**Notes:**

- RNA-seq files must be **gzip-compressed** FASTQ files (`.fastq.gz`).
- PacBio HiFi reads must be provided as a BAM file.
- The pipeline will automatically detect and process multiple RNA-seq samples if present.

---

## 2. Configuration File (`config/config.yaml`)

Edit `config/config.yaml` to match your dataset and analysis requirements.  
Below are the available parameters:

| Parameter               | Description                                                  | Example                                    |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------ |
| `ploidy`                | Ploidy level of the genome (used by GenomeScope2).           | `"2"` for diploid                          |
| `fcs_gx_taxid`          | NCBI Taxonomy ID for FCS-GX screening. [NCBI Taxonomy Tree](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree/) | `"122299"` for *Dioncophyllum thollonii* |
| `busco_lineage_dataset` | BUSCO lineage dataset for genome completeness assessment. [Lineage list](https://busco-data.ezlab.org/v5/data/lineages/) | `"embryophyta_odb12"`                      |
| `dfam_version`          | Version of the Dfam database for RepeatMasker. [Dfam releases](https://www.dfam.org) | `"3.9"`                                    |
| `dfam_lineage_id`       | Dfam lineage ID. See [README.txt](https://www.dfam.org/releases/current/families/FamDB/README.txt). | `"6"` (Viridiplantae)                      |
| `dfam_lineage_name`     | Name of the Dfam lineage to use.                             | `"Viridiplantae"`                          |
| `orthodb_version`       | Version of the OrthoDB database (used by Braker3). [ProtHint instructions](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) | `"12"`                                     |
| `orthodb_lineage`       | OrthoDB lineage dataset to use. [Lineage list](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"Viridiplantae"`                          |
| `orthodb_md5sum`        | MD5 checksum of the OrthoDB database. [Checksums](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/) | `"34c1f027a1a7b10f225b69fbd5500587"`       |
