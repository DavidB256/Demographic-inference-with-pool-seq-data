#!/bin/bash

#SBATCH --mem 9G
#SBATCH --time 00:06:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/data/slurm_errors/slurm-%j.err
#SBATCH -o /scratch/djb3ve/data/slurm_output/slurm-%j.out

# Load conda environment containing moments
source activate myenv

# Handle command line arguments
if [ $# -le 3 ]; then
  echo "Error: Three command line arguments required."
  exit
fi
input_sfs_file=$1
optimization_repeats=$2
killswitch_minutes=$3
ploidy=2
if [ $# -ge 4 ]; then
  ploidy=$4
fi

script="/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/first_pipeline/run_moments_on_pooled_sfs.py"
output_file="/scratch/djb3ve/data/first_models/moments_output.txt"
python $script $input_sfs_file $optimization_repeats $killswitch_minutes $ploidy >> $output_file
