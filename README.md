# genome_assembly_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io) [![Tests](https://github.com/mkrg01/genome_assembly_pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/mkrg01/genome_assembly_pipeline/actions/workflows/ci.yml)

An integrated pipeline for eukaryotic genome assembly and gene annotation. Assembly requires PacBio HiFi reads; ultra-long Oxford Nanopore (ONT) and paired-end Hi-C reads are optional. Gene prediction uses local paired-end RNA-seq or public SRA data.

## Requirements

- [Snakemake ≥ 9.0.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Apptainer (Singularity)](https://apptainer.org/docs/admin/main/installation.html)

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/mkrg01/genome_assembly_pipeline.git
cd genome_assembly_pipeline
```

### 2. Prepare Inputs and Configure Settings

Prepare input files and edit `config/config.yml` following the [configuration guide](config/README.md).

### 3. Run the Workflow

Run from the repository root, adjusting `--cores` to your available resources:

```bash
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 64 all
```

Outputs are written to `results/`. See the [output directory guide](docs/output_directory_structure.md) for details.

## Stepwise Execution

Replace `all` with a target below. Required upstream steps run automatically.

| Target | Description |
| --- | --- |
| `assembly_all` | Hifiasm assembly and QC. |
| `remove_organelle_all` | Organelle removal and QC. |
| `remove_contamination_all` | FCS contamination removal and QC. |
| `purge_dups_all` | Optional haplotig removal and QC. |
| `longstitch_all` | LongStitch correction/scaffolding. |
| `scaffold_all` | YaHS scaffolding and Juicebox contact maps when Hi-C reads are configured. |
| `renamed_all` | Sequence renaming by decreasing length and QC. |
| `softmask_all` | RepeatMasker softmasking. |
| `gene_prediction_all` | Gene prediction, QC, and nuclear release files. |
| `gene_synteny_all` | Gene-based self-synteny and plots. |
| `circos_plot_all` | Gene prediction, QC, release files, and Circos/linear plots. |
| `organelle_annotation_all` | Organelle annotation, circular maps, and release files. |
| `annotation_all` | Nuclear annotation, QC, release files, synteny, and plots. |
| `plots_all` | Circos/linear plots and their required inputs. |
