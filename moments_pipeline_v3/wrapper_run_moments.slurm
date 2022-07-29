#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 1-00:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/moments_pipeline_v3/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/moments_pipeline_v3/output_slurm/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/moments_pipeline_v3/instructions_for_run_moments.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/moments_pipeline_v3/run_moments.py"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
