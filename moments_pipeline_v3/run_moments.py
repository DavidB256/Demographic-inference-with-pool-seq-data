import moments
import numpy as np
import sys
import re
import os
import yaml

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

def two_pop_split(params, ns):
    # params=[nu1, nu2, T, m]
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

def initialize_output_file(output_file):
    with open(output_file, "w") as f:
        column_names = ["sample_size",
                        "mseed",
                        "depth",
                        "pseed",
                        "max_ll"
                        "theta",
                        "n_1_est",
                        "n_2_est",
                        "t_split_est",
                        "mig_rate_est"]
        header_string = ""
        for column_name in column_names:
            header_string += column_name + " \t"
        f.write(header_string + "\n")

if __name__ == "__main__":
    # Handle command line arguments
    sfs_name = sys.argv[1]
    # Use regex to extract numbers from string "sfs_file".
    sfs_name_params = re.findall(r'\d+', sfs_name)

    # Create output file and write its header if it does not already exist
    if not os.path.exists(yd["pipeline_params"]["output_file"]):
        initialize_output_file(yd["pipeline_params"]["output_file"])

    sfs = moments.Spectrum(np.load(yd["pipeline_params"]["data_dir"] + "sfss/" + sfs_name),
                           mask_corners=True, pop_ids=["pop0", "pop1"])

    # Setup for moments
    lower_bound = yd["moments_params"]["lower_bound"]
    upper_bound = yd["moments_params"]["upper_bound"]
    ns = [(int(sfs_name_params[0]) * yd["dem_params"]["ploidy"]) for i in range(2)]

    # Perform moments optimization "optimizaiton_iterations" many times in order
    # to avoid the local extrema problem.
    # Keys are log-likelihoods, values are lists of the form [theta, n_1_est, n_2_est, t_split_est, mig_rate_est]
    moments_output_dict = {}
    for i in range(yd["moments_params"]["num_of_inference_repeats"]):
        params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
        popt = moments.Inference.optimize_log(params, sfs, two_pop_split,
                                              lower_bound=lower_bound,
                                              upper_bound=upper_bound,
                                              maxiter=yd["moments_params"]["max_iter"],
                                              verbose=0)
        model = two_pop_split(popt, ns)
        log_likelihood = moments.Inference.ll_multinom(model, sfs)
        theta = moments.Inference.optimal_sfs_scaling(model, sfs)
        moments_output_dict[log_likelihood] = [theta, popt[0], popt[1], popt[2], popt[3]]

    # Get moments' estimates from iteration with greatest likelihood
    max_log_likelihood = sorted(moments_output_dict)[-1]
    output_list = sfs_name_params + [str(max_log_likelihood)]
    # Convert units from moments units, which depend on theta, mutation rate, and genome length, to conventional units
    # optimal_moments_output = [theta, nu1, nu2, T, m]
    optimal_moments_output = moments_output_dict[max_log_likelihood]
    theta = optimal_moments_output[0]
    conversion_coeff = theta / (4 * yd["dem_params"]["mutation_rate"] * yd["dem_params"]["seq_len"])
    optimal_moments_output[1] *= conversion_coeff
    optimal_moments_output[2] *= conversion_coeff
    optimal_moments_output[3] *= 2 * conversion_coeff
    optimal_moments_output[4] /= 2 * conversion_coeff

    # Export data
    for output in optimal_moments_output:
        output_list.append(str(optimal_moments_output))
    with open(yd["pipeline_params"]["output_file"], "a") as f:
        for output in output_list:
            f.write(output + "\t")
            print(output, end="\t")
        f.write("\n")





