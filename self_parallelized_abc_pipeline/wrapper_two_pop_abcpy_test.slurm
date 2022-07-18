#!/bin/bash

#SBATCH --mem 3G
#SBATCH --time 10:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/self_parallelized_abc_pipeline/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/self_parallelized_abc_pipeline/output_slurm/slurm-%a.err

# Activate environment with moments installed
source activate myenv

script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/self_parallelized_abc_pipeline/two_pop_abcpy_test.py"
options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/self_parallelized_abc_pipeline/threshold_options.txt"

# OPTS = [sfs_name]
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
