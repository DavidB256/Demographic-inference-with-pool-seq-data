import moments
import numpy as np
import sys
import re
import statistics
import os

def two_pop_split(params, ns):
    # params=[nu1, nu2, T, m]
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

def main(demography_model, vcf_file, popinfo_file, vcf_name_params, num_of_optimization_repeats, mutation_rate, genome_length, ploidy):
    output_string = demography_model + "\t"
    for param in vcf_name_params:
        output_string += param + "\t"
    # Placeholder zeros denote that pool-seq noise has not been applied
    output_string += "0\t" * 2
    output_string += str(ploidy) + "\t"

    if demography_model == "two_pop_split":
        ns = [(int(vcf_name_params[0]) * ploidy) for i in range(2)]
        data_dict = moments.Misc.make_data_dict_vcf(vcf_file, popinfo_file)
        sfs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"], projections=ns, polarized=False)

        # Setup for moments
        lower_bound = [1e-3 for i in range(4)]
        upper_bound = [1 for i in range(4)]

        # Perform moments optimization "optimizaiton_iterations" many times in order
        # to avoid the local extrema problem.
        # Keys are log-likelihoods, values are lists of the form [theta, n_1_est, n_2_est, t_split_est, mig_rate_est]
        moments_output_dict = {}
        for i in range(num_of_optimization_repeats):
            params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
            popt = moments.Inference.optimize(params, sfs, two_pop_split,
                                              lower_bound=lower_bound,
                                              upper_bound=upper_bound,
                                              maxiter=10, verbose=100,
                                              multinom=True)
            model = two_pop_split(popt, ns)
            log_likelihood = moments.Inference.ll_multinom(model, sfs)
            theta = moments.Inference.optimal_sfs_scaling(model, sfs)
            moments_output_dict[log_likelihood] = [theta, popt[0], popt[1], popt[2], popt[3]]

        # Get moments' estimates from iteration with greatest likelihood
        max_log_likelihood = sorted(moments_output_dict)[-1]
        output_string += str(max_log_likelihood) + "\t"
        # Convert units from moments units, which depend on theta, mutation rate, and genome length, to conventional units
        optimal_moments_output = moments_output_dict[max_log_likelihood]
        print(optimal_moments_output)
        theta = optimal_moments_output[0]
        conversion_coeff = theta / (4 * mutation_rate * genome_length)
        optimal_moments_output[1] *= conversion_coeff
        optimal_moments_output[2] *= conversion_coeff
        optimal_moments_output[3] *= 2 * conversion_coeff
        optimal_moments_output[4] /= 2 * conversion_coeff
        # Add converted units to "output_string"
        for field in optimal_moments_output:
            output_string += str(field) + "\t"

        return output_string

# Environment variables
output_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/moments_output_from_vcfs.txt"
# Create output file and write its header, if it does not already exist
if not os.path.exists(output_file):
    with open(output_file, "w") as f:
        column_names = ["demography_model",
                        "sample_size_per_pop",
                        "msprime_seed",
                        "poolseq_depth",
                        "poolseq_seed",
                        "ploidy",
                        "max_log_likelihood"
                        "theta",
                        "n_1_est",
                        "n_2_est",
                        "t_split_est",
                        "mig_rate_est"]
        header_string = ""
        for column_name in column_names:
            header_string += column_name + "\t"
        f.write(header_string + "\n")
input_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"

# Handle command line arguments, which come from "pipeline_instructions.txt"
vcf_name = sys.argv[1]
popinfo_name = sys.argv[2]
num_of_optimization_repeats = int(sys.argv[3])
mutation_rate = float(sys.argv[4])
genome_length = int(sys.argv[5])
ploidy = int(sys.argv[6])
# Use regex to extract numbers from string "sfs_file".
# sfs_name_params = [n, msprime_seed, depth, poolseq_seed]
vcf_name_params = re.findall(r'\d+', vcf_name)
# Determine demography model type
if vcf_name.startswith("two_pop_split"):
    demography_model = "two_pop_split"
else:
    print("Error: Unrecognized demography model.")
    sys.exit()

output_string = main(demography_model, input_dir + "vcfs/" + vcf_name, input_dir + "popinfos/" + popinfo_name,
                     vcf_name_params, num_of_optimization_repeats, mutation_rate, genome_length, ploidy)
with open(output_file, "a") as f:
    f.write(output_string + "\n")
