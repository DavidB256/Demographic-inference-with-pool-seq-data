#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 1:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/pmcabc_pipeline/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/pmcabc_pipeline/output_slurm/slurm-%a.err

# Activate environment with moments installed
source activate myenv

script_file=/scratch/djb3ve/Demographic-inference-with-pool-seq-data/pmcabc_pipeline/two_pop_abcpy_test.py
options_file=/scratch/djb3ve/Demographic-inference-with-pool-seq-data/pmcabc_pipeline/threshold_options.txt

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
