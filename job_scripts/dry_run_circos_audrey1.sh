#!/bin/bash
# SLURM on audrey1, a remote server at Fukushima Lab

#SBATCH -J dry_run_genome_assembly_pipeline
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=16G # Memory per CPU (max 480 GB on this cluster)
#SBATCH -t 2976:00:00 # Maximum time in d-hh:mm:ss format.
#SBATCH --output=dry_run_genome_assembly_pipeline_%A_%a.out
#SBATCH --error=dry_run_genome_assembly_pipeline_%A_%a.err
#SBATCH -p debug # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

echo "`date`: Starting"

# Activate Conda environment
source ~/miniforge3/etc/profile.d/conda.sh # Modify the path according to your installation
conda activate snakemake
echo "Snakemake version: `snakemake --version`"

# Dry run Snakemake workflow
snakemake -n --sdm apptainer --singularity-args "--bind $(pwd)" --cores 1 -s workflow/Snakefile.circos all

echo "`date`: Ending"
