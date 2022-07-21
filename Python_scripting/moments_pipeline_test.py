import Optimize_Functions
import moments
import numpy as np
import sys

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/alt_2island_3mig_model.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
ns = [20, 20]
prefix = sys.argv[1]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def null_split(params, ns):
    nu1, nu2, T, mig_rate = params
    return moments.Demographics2D.split_mig([nu1, nu2, T, mig_rate], ns,
                                            pop_ids=["pop0", "pop1"])

Optimize_Functions.Optimize_Routine(fs=fs, outfile=prefix,
                                    model_name="null_split",
                                    func=null_split,
                                    rounds=3, param_number=4,
                                    fs_folded=True,
                                    maxiters=[100 for i in range(3)],
                                    param_labels=["nu1", "nu2", "T_split", "mig_rate"])
