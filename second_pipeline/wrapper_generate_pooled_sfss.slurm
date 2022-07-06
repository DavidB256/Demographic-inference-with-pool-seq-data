#!/bin/bash

#SBATCH --mem 128G
#SBATCH --time 10:00:00
#SBATCH --partition largemem
#SBATCH --account biol4585j-yey2sn
#SBATCH -c 1
#SBATCH -e /scratch/djb3ve/data/r_slurm_errors/slurm-%a.err
#SBATCH --output /scratch/djb3ve/data/r_slurm_output/slurm-%a.out

module load gcc/7.1.0 openmpi/3.1.4 intel/18.0 intelmpi/18.0 R/4.1.1

options_file="/scratch/djb3ve/data/second_pipeline/pipeline_instructions.txt"
script_file="/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/second_pipeline/generate_pooled_sfss.R"

# OPTS = [vcf_file, popinfo_eval_string, haploid_counts_eval_string, poolseq_depth, popinfo_file]
# The last argument is not considered by this R script. It only exists for use of moments_output
# later down the pipeline.
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${options_file})

Rscript $script_file $OPTS
