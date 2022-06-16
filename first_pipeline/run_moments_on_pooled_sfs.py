import moments
import numpy as np
import sys
import os
import re
import statistics

def two_pop_split(params, ns):
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

# Returns numpy array containing SFS imported from serialized file "file_name"
def import_sfs_from_serialized(file_name):
    with open(file_name, "r") as f:
        # Read header of serialized file, which contains dimensions of SFS
        dims = []
        line = f.readline()
        while line[:-1].isnumeric():
            dims.append(int(line[:-1]))
            line = f.readline()

        # Read elements of SFS from body of serialized file
        values = []
        line = f.readline()
        while line:
            values.append(int(line[:-1]))
            line = f.readline()

    # Place values into numpy array
    sfs = np.transpose(np.array(values).reshape(dims))
    # This hardcoding for the two_pop_split model needs to be replaced
    return moments.Spectrum(sfs, mask_corners=True, pop_ids=["pop0", "pop1"])

def run_moments_on_two_pop_split(sfs, sfs_params, iterations):
    lower_bound = [1e-3 for i in range(4)]
    upper_bound = [1e-2 for i in range(4)]

    ns = [int(sfs_params[0]) * ploidy for i in range(2)]

    lls = []
    for i in range(iterations):
        print(i)
        params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
        print(sfs)
        popt = moments.Inference.optimize(params, sfs, two_pop_split,
                                              lower_bound=lower_bound,
                                              upper_bound=upper_bound,
                                              maxiter=10, verbose=0,
                                              multinom=True)
        model = moments.Demographics2D.split_mig(popt, ns)
        lls.append(moments.Inference.ll_multinom(model, sfs))

    col_values = ["two_pop_split"] + sfs_params + [str(statistics.median(lls)), str(max(lls))]
    return_string = ""
    for value in col_values:
        return_string += value + "\t"
    print(return_string)
    return return_string + "\n"


# Global variables
input_dir = "/scratch/djb3ve/data/first_models/serialized_pooled_sfss/"
output_file = "/scratch/djb3ve/data/first_models/moments_output/pooled_sfss.txt"
iterations = 5
ploidy = 2

with open(output_file, "w") as f:
    f.write("model\tn\tmsprime_seed\tdepth\tpoolseq_seed\tmedian_ll\tmax_ll")
    for sfs_file in os.listdir(input_dir):
        if sfs_file != "two_pop_split_demography_n50_seed9_depth100_pooled_sfs_serialized_seed1.txt":
            continue
        print(sfs_file)
        if sfs_file.startswith("two_pop_split_demography"):
            sfs = import_sfs_from_serialized(input_dir + sfs_file)
            # sfs_params = [n, msprime_seed, depth, poolseq_seed]
            sfs_params = re.findall(r'\d+', sfs_file)
            f.write(run_moments_on_two_pop_split(sfs, sfs_params, iterations))
