#! /bin/bash

# AGE on SHIROKANE, the supercomputer at the Human Genome Center
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4 # Number of CPUs
#$ -l s_vmem=16G # Memory per CPU
#$ -l mjob # https://gc.hgc.jp/uge/uge-resource/

echo "`date`: Starting"

# Activate Apptainer
module use /usr/local/package/modulefiles
module load apptainer
echo "Apptainer version: `apptainer --version`"

# Activate Conda environment
source ~/miniforge3/etc/profile.d/conda.sh # Modify the path according to your installation
conda activate snakemake
echo "Snakemake version: `snakemake --version`"

# Dry run Snakemake workflow
snakemake -n --sdm apptainer --singularity-args "--bind $(pwd)" --cores 1

echo "`date`: Ending"
