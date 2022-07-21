#!/bin/bash

#SBATCH --mem 1G
#SBATCH --time 1:00:00
#SBATCH --partition standard
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/slurm_errors/slurm-%a.err
#SBATCH --output /scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/slurm_output/slurm-%a.out

options_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/pipeline_instructions.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-pool-seq-data/second_pipeline/generate_pooled_sfss.py"
num_of_replicates=3

# OPTS = [vcf_name, popinfo_name, haploid_counts_eval_string, poolseq_depth, rounding_method, ploidy]
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})
python $script_file $OPTS $num_of_replicates
