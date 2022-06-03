import Optimize_Functions
import moments
import numpy as np
import sys

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/alt_2island_1mig_model.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
ns = [20, 20]
prefix = sys.argv[1]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def two_island_admixture(params, ns):
    nu1, nu2, T = params
    mig_rate = 0
    return moments.Demographics2D.split_mig([nu1, nu2, T, mig_rate], ns,
                                            pop_ids=["pop0", "pop1"])

    return fs

Optimize_Functions.Optimize_Routine(fs, prefix, "two_island_admixture",
                                    two_island_admixture, 3, 3, fs_folded=True)
