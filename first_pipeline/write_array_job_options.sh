#!/bin/bash

file_dir="/scratch/djb3ve/data/first_models/serialized_pooled_sfss/"
output_file="/scratch/djb3ve/data/first_models/array_job_options.txt"

# Delete output file if it already exists
if [ -f $output_file ]; then
  rm $output_file
fi

for file in $file_dir*; do
  echo -e "${file}\t5\t2" >> $output_file
done
