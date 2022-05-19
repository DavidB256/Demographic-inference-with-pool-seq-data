#!/bin/bash

rr_exps="6 7 8 9 10"
mu_exps="6 5 4 3 2"
out="combined_pi_output.txt"

if [ -f $out ]; then
	rm $out
fi

for rr_exp in $rr_exps; do
	for mu_exp in $mu_exps; do
		cat dir_pi_out/${rr_exp}rr_${mu_exp}mu_out.vcf.windowed.pi | tail -n +2 | while read line; do echo -e "${line}\t${rr_exp}\t${mu_exp}"; done >> $out
	done
done
