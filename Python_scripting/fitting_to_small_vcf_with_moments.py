import moments
import dadi
import numpy as np
import sys

iterations = int(sys.argv[1])

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/small_vcf.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_for_small_vcf.txt"
output = "/scratch/djb3ve/Demographic-inference-with-Pool-seq-data/moments_fitting_small.txt"
ns = [5]

print("Setup complete.")

data_dict = dadi.Misc.make_data_dict_vcf(vcf, popinfo)

print("VCF imported.")

fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"],
                                     projections=ns, polarized=False)

print("VCF converted to SFS.")

def isolated_island(params, ns):
    nu1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    fs.integrate([nu1], 100000)
    return fs

print("Model function defined.")

lower_bound = [1e-3]
upper_bound = [1000]

out_f = open(output, "w")

for i in range(iterations):
    params = [np.random.uniform(lower_bound[j], upper_bound[j])
              for j in range(1)]
    popt = moments.Inference.optimize_log(params, fs, isolated_island,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound)
    model = isolated_island(popt, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)

    output_string = "%f\t%f\t%f\n" % (params[0], popt[0], ll_model)
    print(output_string, end="")
    out_f.write(output_string)

out_f.close()
