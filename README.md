# genome_assembly_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io) [![Tests](https://github.com/mkrg01/genome_assembly_pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/mkrg01/genome_assembly_pipeline/actions/workflows/ci.yml)

This is an integrated pipeline for eukaryotic genome assembly and gene annotation. PacBio HiFi reads are required for the full workflow. Paired-end RNA-seq reads are required for gene prediction. Ultra-long Oxford Nanopore (ONT) reads and paired-end Hi-C reads are optionally supported.

If you already have an assembly FASTA from another workflow, use `workflow/Snakefile.annotation` to run the downstream annotation path only: RepeatModeler/RepeatMasker softmasking, BRAKER3 gene prediction, protein/transcript QC, formatting, and Circos/linear plots. In that mode, set `external_assembly` in `config/config.yml` and keep paired-end RNA-seq reads in `raw_data/`.

See [this page](docs/output_directory_structure.md) for details on the expected outputs.

## Requirements

Before running the workflow, make sure the following software is installed:

- [Snakemake ≥ 9.0.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Apptainer (Singularity)](https://apptainer.org/docs/admin/main/installation.html)

## Getting Started

Follow the steps below to set up and run the workflow:

### 1. Clone the Repository

Clone this repository to your local machine:

```
git clone https://github.com/mkrg01/genome_assembly_pipeline.git
cd genome_assembly_pipeline
```

### 2. Prepare Input Files and Configure Settings

See [`config/README.md`](config/README.md) for details on preparing input files and adjusting configuration parameters.
The configuration guide follows the same section order as `config/config.yml`, from runtime settings through assembly, contamination screening, YaHS scaffolding, QC, visualization, and organelle annotation options.

### 3. Execute the Workflow

Run the workflow from the repository root directory. Replace `/path/to/repo` with the actual path to your local repository:

```
cd /path/to/repo
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 64 all
```

> [!TIP]
> You can run the pipeline in a stepwise manner. Replace `all` with one of the commands below.
> 
> 1. `assembly_all`: Runs rules up to the generation of the Hifiasm assembly and its associated metrics.
> 2. `remove_organelle_all`: Runs rules up to the organelle removal step and its associated metrics.
> 3. `remove_contamination_all`: Runs rules up to the contamination removal step by FCS and its associated metrics.
> 4. `scaffold_all`: Runs rules up to YaHS Hi-C scaffolding and Juicebox-ready contact maps when Hi-C reads are configured. Without Hi-C reads, this is effectively the same as `remove_contamination_all`.
> 5. `softmask_all`: Runs rules up to softmasking by RepeatMasker.
> 6. `gene_prediction_all`: Runs rules up to gene prediction and its associated metrics.
> 7. `circos_plot_all`: Runs rules up to the Circos plot for the main genome analysis path.
> 8. `organelle_annotation_all`: Annotates Oatk-assembled organelle genomes, draws pyCirclize and gbdraw circular maps, and stages organelle genome and annotation files for the release package.
> 
> You do not need to start from step 1 — for example, if you run `remove_contamination_all` first, the rules related to `assembly_all` and `remove_organelle_all` will be executed automatically.

For organelle annotation curation, a typical sequence is:

```bash
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 64 organelle_annotation_all
# Manually curate results/organelle_annotation/.../*.gbk if necessary.
snakemake -n --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 4 organelle_annotation_all --rerun-triggers mtime
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 4 organelle_annotation_all --rerun-triggers mtime
```

After manual GenBank curation, the dry-run should schedule only downstream steps such as organelle map drawing and release staging.
If GenBank-producing rules are listed, check the upstream file timestamps before running the command without `-n`.

> [!NOTE]
> Adjust the `--cores` value based on your available computational resources.

> [!NOTE]
> All rules except those with FCS wrapper scripts (`fcs.py`, `run_fcsadaptor.sh`) run in containers. These wrapper scripts internally call the main FCS functions, which are executed inside containers.

The output will be generated in the [`results` directory](docs/output_directory_structure.md).

## Annotating an External Assembly

To annotate an existing genome assembly, set `external_assembly` in `config/config.yml` and run:

```
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 64 -s workflow/Snakefile.annotation all
```
