#!/bin/bash

#SBATCH --mem 2G
#SBATCH --time 1:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/errors_slurm/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/output_slurm/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/instructions/instructions_for_generate_sfss.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/third_pipeline/generate_sfss.py"

# OPTS = [vcf_name, popinfo_name, haploid_counts_eval_string, poolseq_depth, rounding_method, ploidy]
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS
