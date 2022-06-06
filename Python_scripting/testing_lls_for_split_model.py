import moments
import numpy as np
import sys


# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/msprime_null_split_small.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
ns = [20, 20]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def two_island_admixture(params, ns):
    # params = [nu1, nu2, T_split, mig_rate]
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

lower_bound = [1e-2 for i in range(4)]
upper_bound = [1e2 for i in range(4)]

def fit_params(starting_params):
    # popt = moments.Inference.optimize_log(starting_params, fs, two_island_admixture,
    #                                       lower_bound=lower_bound,
    #                                       upper_bound=upper_bound,
    #                                       maxiter=100, verbose=10)
    # model = two_island_admixture(popt, ns)
    model = two_island_admixture(starting_params, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)
    print("""
          Starting parameters: {0}\n
          LL: {2}\n
          theta: {3}
          """.format(starting_params, ll_model, theta))

params = [1 for i in range(4)]
fit_params(params)
params = [30 / 4e-3, 20 / 4e-3, 1000 / 2e-3, 3e-5 / 2e-3]
fit_params(params)
