#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 1-00:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline_v1/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline_v1/output_slurm/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline_v1/instructions_for_run_abc.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline_v1/run_abc.py"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS



