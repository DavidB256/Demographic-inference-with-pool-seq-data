#!/bin/bash

out="combined_models.pi"
source="/scratch/djb3ve/msprime/may16/dir_model_pi_out/"

if [ -f $out ]; then
	rm $out
fi

cat ${source}control.vcf.windowed.pi | tail -n +2 | while read line; do echo -e "${line}\tcontrol"; done >> $out
cat ${source}continent_islands.vcf.windowed.pi | tail -n +2 | while read line; do echo -e "${line}\tcontinent"; done >> $out
cat ${source}equal_islands.vcf.windowed.pi | tail -n +2 | while read line; do echo -e "${line}\tislands"; done >> $out
cat ${source}stepping_stone.vcf.windowed.pi | tail -n +2 | while read line; do echo -e "${line}\tstepping"; done >> $out
