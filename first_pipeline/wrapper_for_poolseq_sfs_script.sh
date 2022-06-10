#!/bin/bash

# Iterates through VCF files created by generate_demography_vcfs_for_initial_models
# and outputs serialized SFS outputs from poolseq_sfs_script.R into a subdirectory of
# /scratch/djb3ve/data/

module purge
module load gcc/7.1.0 openmpi/3.1.4 intel/18.0 intelmpi/18.0 R/4.1.1

# Handle command line arguments
# First required argument is input directory "input_dir"
# Second required argument is output directory "output_dir"
# Third optional argument is random seed "random_seed", which defaults to 1
if [ $# -ge 2 ]; then
  input_dir=$1
  output_dir=$2
  if [ $# -ge 3 ]; then
    random_seed=$3
  else
    random_seed=1
  fi
else
  echo "Error: At least two command line arguments are required."
  exit
fi

for vcf in $input_dir*; do
  Rscript /scratch/djb3ve/first_pipeline/poolseq_sfs_script.R $vcf $output_dir $random_seed
done
