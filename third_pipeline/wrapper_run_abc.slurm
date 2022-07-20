#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 12:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/output_slurm/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/instructions/instructions_for_run_abc.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/run_abc.py"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS



