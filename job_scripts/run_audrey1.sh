#!/bin/bash
# SLURM on audrey1, a remote server at Fukushima Lab

#SBATCH -J genome_assembly_pipeline
#SBATCH -c 128 # Number of CPUs
#SBATCH --mem-per-cpu=3.75G # Memory per CPU (max 480 GB on this cluster)
#SBATCH -t 2976:00:00 # Maximum time in d-hh:mm:ss format.
#SBATCH --output=genome_assembly_pipeline_%A_%a.out
#SBATCH --error=genome_assembly_pipeline_%A_%a.err
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

# Run Snakemake workflow
snakemake --sdm conda apptainer --singularity-args "--bind $(pwd)" --cores "$SLURM_CPUS_PER_TASK"

echo "`date`: Ending"
