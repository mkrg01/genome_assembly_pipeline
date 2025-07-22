
# repository_name

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

The workflow is managed using [Snakemake](https://snakemake.github.io/).

## Getting Started

Follow these steps to set up and run the workflow:

### 1. Clone the Repository

Clone this repository to your local machine:

```
git clone https://github.com/mkrg01/repository_name.git

cd repository_name
```

### 2. Install Snakemake

Ensure you have [Snakemake version 8 or higher installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

### 3. Execute the Workflow

Run the workflow from the repository’s root directory. Replace `/path/to/repo` with the actual path to the repository on your system:

```
cd /path/to/repo

snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores 8
```

**Note:** Adjust the `--cores` parameter based on your available computational resources. 

The output files will be generated in the `results` folder.

## License

This repository is licensed under the [MIT license](LICENSE).