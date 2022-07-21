#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 10:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/b_slurm_errors/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/b_slurm_output/slurm-%a.out

# Activate environment with moments installed
source activate myenv

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/sfs_name_list.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/second_pipeline/run_moments_on_pooled_sfs.py"
num_of_optimization_repeats=3
# 2.9e-9
mutation_rate=0.0000000029
# 1e6
genome_length=1000000
ploidy=2

# OPTS = [sfs_name]
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS $num_of_optimization_repeats $mutation_rate $genome_length $ploidy
