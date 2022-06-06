import moments
import numpy as np
import sys

iterations = int(sys.argv[1])

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/msprime_null_split_small.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
output = "/scratch/djb3ve/data/moments_fitting_split_demography.txt"
ns = [20, 20]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def two_island_admixture(params, ns):
    # params = [nu1, nu2, T_split, mig_rate]
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

lower_bound = [1e-3 for i in range(4)]
upper_bound = [1e2 for i in range(4)]

out_f = open(output, "w")
header_string = "nu1_initial\t" + \
                "nu2_initial\t" + \
                "T_initial\t" + \
                "m_initial\t" + \
                "nu1_optimized\t" + \
                "nu_2_optimized\t" + \
                "T_optimized\t" + \
                "m_optimized\t" + \
                "theta\t" + \
                "log-likelihood\n"
out_f.write(header_string)
print(header_string, end="")

for i in range(iterations):
    print("here1")
    params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
    print("here2")
    popt = moments.Inference.optimize_log(params, fs, two_island_admixture,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound,
                                          maxiter=100, verbose=1)
    print("here3")
    model = two_island_admixture(popt, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)
    output_string = "%f\t" % params[0] + \
                    "%f\t" % params[1] + \
                    "%f\t" % params[2] + \
                    "%f\t" % params[3] + \
                    "%f\t" % popt[0] + \
                    "%f\t" % popt[1] + \
                    "%f\t" % popt[2] + \
                    "%f\t" % popt[3] + \
                    "%f\t" % theta + \
                    "%f\t" % ll_model + "\n"
    out_f.write(output_string)
    print(output_string, end="")

out_f.close()
