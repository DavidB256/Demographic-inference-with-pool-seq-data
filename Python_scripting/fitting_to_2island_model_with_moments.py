import moments
import dadi
import numpy as np
import sys

iterations = int(sys.argv[1])

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/2island_1mig_model.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
output = "/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/moments_fitting_2islands.txt"
ns = [1, 1]

print("Setup complete.")

data_dict = dadi.Misc.make_data_dict_vcf(vcf, popinfo)

print("VCF imported.")

fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=ns, polarized=False)

print("VCF converted to SFS.")

def two_island_admixture(params, ns):
    nu1, nu2, T = params
    mig_rate = 0.1

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, mig_rate], [mig_rate, 0]]))

    return fs

print("Model function defined.")

lower_bound = [90, 90, 100000]
upper_bound = [110, 110, 100000]

out_f = open(output, "w")

for i in range(iterations):
    params = [np.random.uniform(lower_bound[j], upper_bound[j])
              for j in range(3)]
    popt = moments.Inference.optimize_log(params, fs, two_island_admixture,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound)
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
