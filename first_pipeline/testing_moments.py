import moments
import numpy as np
import sys

iterations = int(sys.argv[1])

vcf = "/scratch/djb3ve/data/first_models/two_pop_split_demography_n10_seed1.vcf"
popinfo = "/scratch/djb3ve/data/first_models/two_pop_split_demography_n10_popinfo.txt"
ns = [10, 10]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"], projections=ns, polarized=False)

lower_bound = [1e-3 for i in range(4)]
upper_bound = [1e-1 for i in range(4)]

column_names_list = ["nu1_init", "nu2_init", "t_init", "mig_init",
                     "nu1_opt", "nu2_opt", "t_opt", "mig_opt",
                     "ll", "theta"]
header = ""
for name in column_names_list:
    header += name + "\t"
print(header)

for i in range(iterations):
    params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
    popt = moments.Inference.optimize_log(params, fs, two_pop_split_model,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound,
                                          maxiter=100, verbose=0,
                                          multinom=True)
    model = two_pop_split_model(popt, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)

    print("%f\t" * len(column_names_list) %
          (params[0], params[1], params[2], params[3],
           popt[0], popt[1], popt[2], popt[3],
           ll_model, theta))
