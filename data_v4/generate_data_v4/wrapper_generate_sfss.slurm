#!/bin/bash

# Adapted for data_v4.

#SBATCH --mem 3G
#SBATCH --time 1:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data_v4/generate_data_v4/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data_v4/generate_data_v4/output_slurm/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data_v4/generate_data_v4/instructions_for_generate_sfss.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data_v4/generate_data_v4/generate_sfss.py"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
