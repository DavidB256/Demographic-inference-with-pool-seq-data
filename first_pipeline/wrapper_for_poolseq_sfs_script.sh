#!/bin/bash

# Iterates through VCF files created by generate_demography_vcfs_for_initial_models
# and outputs serialized SFS outputs from poolseq_sfs_script.R into a subdirectory of
# /scratch/djb3ve/data/

module purge
module load gcc/7.1.0 openmpi/3.1.4 intel/18.0 intelmpi/18.0 R/4.1.1
