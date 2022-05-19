#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 24:00:00
#SBATCH -o /home/ksl2za/ABC/errors_outs/outs/sc_out.out
#SBATCH -e /home/ksl2za/ABC/errors_outs/errors/sc_error.error
#SBATCH -p standard
#SBATCH -A lfg_lab

echo "began at"  `date`

#Load conda module
module purge
module load anaconda/2020.11-py3.8

#Activate dadi kernel
source activate dadi_env

ABC_sc=/home/ksl2za/ABC/ABC_sc.py

#run 100,000 iterations of the sc model
cd ~/ABC
python $ABC_sc \
100000

#De-Activate kernel
conda deactivate

#Print the time
echo "ended at"  `date`