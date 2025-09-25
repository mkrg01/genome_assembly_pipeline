# genome_assembly_pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)

This is an integrated pipeline for eukaryotic genome assembly and gene annotation.
It currently supports PacBio HiFi reads and RNA-seq reads, both of which are required as inputs.
See [this page](https://github.com/mkrg01/genome_assembly_pipeline/wiki/Directory-structure-in-results) for details on the expected outputs.

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

See [`config/README.md`](https://github.com/mkrg01/genome_assembly_pipeline/blob/main/config/README.md) for details on preparing input files and adjusting configuration parameters.

### 3. Execute the Workflow

Run the workflow from the repository root directory. Replace `/path/to/repo` with the actual path to your local repository:

```
cd /path/to/repo

snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 48
```
> [!NOTE]
> Adjust the `--cores` value based on your available computational resources.

> [!NOTE]
> All rules will be executed in Singularity containers.

The output will be generated in the [`results` directory](https://github.com/mkrg01/genome_assembly_pipeline/wiki/Directory-structure-in-results).
