import numpy as np
import moments

input_file = "/scratch/djb3ve/data/second_pipeline/vcfs/two_pop_split_demography_n10_seed1.vcf"
popinfo = "/scratch/djb3ve/data/second_pipeline/popinfos/two_pop_split_demography_n10_popinfo.txt"
ns = [20, 20]

data_dict = moments.Misc.make_data_dict_vcf(input_file, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"], projections=ns, polarized=False)

print(fs)
