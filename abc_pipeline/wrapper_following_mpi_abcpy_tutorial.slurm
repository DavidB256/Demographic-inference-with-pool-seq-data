#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --time 1:00
#SBATCH --partition parallel
#SBATCH -A biol4585j-yey2sn
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/error.txt
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/output.txt

# Activate environment with moments installed
source activate myenv
module load gcc openmpi
srun python /scratch/djb3ve/Demographic-inference-with-pool-seq-data/abc_pipeline/following_mpi_abcpy_tutorial.py
