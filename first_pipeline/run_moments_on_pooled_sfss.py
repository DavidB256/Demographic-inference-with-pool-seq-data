import moments
import numpy as np
import sys

# Global variables: working directory and pipeline instructions file
working_dir = "/scratch/djb3ve/data/first_models/"
instructions = working_dir + "pipeline_instructions.txt"

# Returns numpy array containing SFS imported from serialized file "file_name"
def import_sfs_from_serialized(file_name):
    with open(working_dir + file_name, "r") as f:
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
    return sfs

def run_moments_on_two_pop_split(serialized_file, ns):
    # Perform inference with moments
    fs = import_sfs_from_serialized()
    lower_bound = [1e-3 for i in range(4)]
    upper_bound = [1 for i in range(4)]

    column_names_list = ["nu1_init", "nu2_init", "t_init", "mig_init",
                         "nu1_opt", "nu2_opt", "t_opt", "mig_opt",
                         "ll", "theta"]
    header = ""
    for name in column_names_list:
        header += name + "\t"

    for i in range(iterations):
    params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
    popt = moments.Inference.optimize_log(params, fs, moments.Demographics2D.split_mig,
                                          lower_bound=lower_bound,
                                          upper_bound=upper_bound,
                                          maxiter=200, verbose=0,
                                          multinom=True)
    model = moments.Demographics2D.split_mig(popt, ns)
    ll_model = moments.Inference.ll_multinom(model, fs)
    theta = moments.Inference.optimal_sfs_scaling(model, fs)

    print("%f\t" * len(column_names_list) %
          (params[0], params[1], params[2], params[3],
           popt[0], popt[1], popt[2], popt[3],
           ll_model, theta))

with open(instructions, "r") as f:
    # Read through pipeline instructions file
    line = f.readline()
    while line:
        # Skip commented-out lines, like the header
        if line[0] == "#":
            continue

        # Extract relevant data from line of pipeline instructions, which also contains
        # information only needed by "poolseq_sfs_script.R"
        words = line.split()
        vcf = words[0]
        ns = words[1].split("=")[1][2:-2].split(",")
        ns = [int(n) for n in ns]
        depth = words[3]
        popinfo = words[4]

        if vcf.startswith("two_pop_split_demography"):
            run_moments_on_two_pop_split(vcf, ns, depth, popinfo)

        # Read in next line before looping
        line = f.readline()
