#!/bin/bash

input_dir="/scratch/djb3ve/data/first_models/serialized_pooled_sfss/"
output_file="/scratch/djb3ve/data/first_models/moments_output.txt"
echo -e "model\tn\tmsprime_seed\tdepth\tpoolseq_seed\tmedian_ll\tmax_ll" > $output_file

for sfs_file in $input_dir*; do
  echo $sfs_file
  sbatch /scratch/djb3ve/Demographic-inference-with-Pool-seq-data/first_pipeline/wrapper_run_moments_on_pooled_sfs.slurm $sfs_file 10 5 2
done
