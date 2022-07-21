import os
import re

working_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"
output_file = working_dir + "moments_vcf_instructions.txt"

with open(output_file, "w") as f:
    for vcf_name in os.listdir(working_dir + "vcfs/"):
        n = re.findall(r'\d+', vcf_name)[0]
        popinfo_name = f"two_pop_split_demography_n{n}_popinfo.txt"
        f.write("%s\t%s\n" % (vcf_name, popinfo_name))
