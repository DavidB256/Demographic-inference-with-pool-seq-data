#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time 1:00
#SBATCH --partition parallel
#SBATCH -A biol4585j-yey2sn
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/error.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/output.out

# Activate environment with moments installed
source activate myenv
module load gcc openmpi

#options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/thresholds.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/mpi_two_pop_abcpy_test.py"
#OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
srun python $script_file 0.1
