import os
import re

# Writes pipeline_instructions.txt for second_pipeline, which only features
# two_pop_split demography models.

poolseq_depths = [10, 40, 70, 100]
rounding_method = "counts"
ploidy = 2
working_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"
output_file = "pipeline_instructions.txt"

with open(working_dir + output_file, "w") as f:
    for vcf_name in os.listdir(working_dir + "vcfs/"):
        for poolseq_depth in poolseq_depths:
            n = re.findall(r'\d+', vcf_name)[0]
            popinfo_name = f"two_pop_split_demography_n{n}_popinfo.txt"
            haploid_counts_string = f"[{ploidy * int(n)},{ploidy * int(n)}]"

            f.write("%s\t%s\t%s\t%d\t%s\t%s\n" %
                    (vcf_name,
                    popinfo_name,
                    haploid_counts_string,
                    poolseq_depth,
                    rounding_method,
                    ploidy))
