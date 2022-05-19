#!/bin/bash

vcf_pool=/scratch/djb3ve/msprime/may16/dir_vcf_out/

for file in $vcf_pool*; do
	name_id=$(echo $file | sed 's|/scratch/djb3ve/msprime/may16/dir_vcf_out/||g')
	
	vcftools \
	--vcf $file \
	--window-pi 100000 \
	--window-pi-step 100000 \
	--out dir_pi_out/$name_id
done
