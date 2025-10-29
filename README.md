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
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 64 all
```

> [!TIP]
> You can run the pipeline in a stepwise manner. Replace `all` with one of the command below.
> 
> 1. `assembly_all`: Runs rules up to the generation of the Hifiasm assembly and its associated metrics.
> 2. `remove_organelle_all`: Runs rules up to the organelle removal step and its associated metrics.
> 3. `remove_contamination_all`: Runs rules up to the contamination removal step by FCS and its associated metrics.
> 4. `softmask_all`: Runs rules up to softmasking by RepeatMasker.
> 5. `gene_prediction_all`: Runs rules up to gene prediction and related metrics (equivalent to `all`).
> 
> You do not need to start from step 1 — for example, if you run `remove_contamination_all` first, the rules related to `assembly_all` and `remove_organelle_all` will be executed automatically.

> [!NOTE]
> Adjust the `--cores` value based on your available computational resources.

> [!NOTE]
> All rules except those with FCS wrapper scripts (`fcs.py`, `run_fcsadaptor.sh`) run in containers. These wrapper scripts internally call the main FCS functions, which are executed inside containers.

The output will be generated in the [`results` directory](https://github.com/mkrg01/genome_assembly_pipeline/wiki/Directory-structure-in-results).
