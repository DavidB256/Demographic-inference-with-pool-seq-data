#!/bin/bash

# Load most recent version of R on Rivanna
module purge
module load gcc/7.1.0 openmpi/3.1.4 intel/18.0 intelmpi/18.0 R/4.1.1

# Handle command line arguments
if [ $# -eq 0 ]; then
  echo "Error: One command line argument is required."
  exit
fi
instructions=$1

script="/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/first_pipeline/poolseq_sfs_script.R"

while read -r line; do
  # Skip lines that are commented out
  if [ ${line:0:1} == "#" ]; then
    continue
  fi
  Rscript $script $line
done < $instructions
