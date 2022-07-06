#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 10:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/data/slurm_errors/slurm-%a.err
#SBATCH --output /scratch/djb3ve/data/slurm_output/slurm-%a.out

# Load conda environment containing moments
source activate myenv

options_file="/scratch/djb3ve/data/first_models/array_job_options.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/first_pipeline/run_moments_on_pooled_sfs.py"
output_file="/scratch/djb3ve/data/first_models/moments_output.txt"

# OPTS = [input_sfs_file, optimization_repeats, ploidy]
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})

python $script_file $OPTS >> $output_file
