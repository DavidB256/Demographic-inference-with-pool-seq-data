import moments
import numpy as np
import sys

iterations = int(sys.argv[1])

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/2island_1mig_model.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
output = "/scratch/djb3ve/data/moments_fitting_2islands.txt"
ns = [20, 20]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

def two_island_admixture(params, ns):
    nu1, nu2, T = params
    mig_rate = 0.1

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, mig_rate], [mig_rate, 0]]))

    return fs

lower_bound = [1e-3, 1e-3, 1e-3]
upper_bound = [200, 200, 1000]

out_f = open(output, "w")
header_string = "nu1_initial\tnu2_initial\tT_initial\tnu1_optimized\tnu_2_optimized\tT_optimized\ttheta\t\tlog-likelihood\n"
print(header_string, end="")
out_f.write(header_string)

for i in range(iterations):
    params = [np.random.uniform(lower_bound[j], upper_bound[j])
              for j in range(3)]
    popt = moments.Inference.optimize_log(params, fs, two_island_admixture,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound,
                                          maxiter=100)
    model = two_island_admixture(popt, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)
    popt = [i * theta for i in popt]
    output_string = "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % \
                    (params[0], params[1], params[2],
                     popt[0], popt[1], popt[2], theta, ll_model)
    print(output_string, end="")
    out_f.write(output_string)

out_f.close()
