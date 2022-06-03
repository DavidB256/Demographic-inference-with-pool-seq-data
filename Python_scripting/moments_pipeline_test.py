import Optimize_Functions
import moments
import numpy as np
import sys

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/msprime_null_split.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_2_samples_of_50.txt"
ns = [100, 100]
prefix = sys.argv[1]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def null_split(params, ns):
    nu1, nu2, T, mig_rate = params
    return moments.Demographics2D.split_mig([nu1, nu2, T, mig_rate], ns,
                                            pop_ids=["pop0", "pop1"])
    return fs

Optimize_Functions.Optimize_Routine(fs, prefix, "null_split",
                                    null_split, 3, 4, fs_folded=True)
