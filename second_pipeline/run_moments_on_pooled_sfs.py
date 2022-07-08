import moments
import numpy as np

# -----
# COPIED FROM SAME VERSION OF THIS SCRIPT FROM FIRST_PIPELINE
# -----

def two_pop_split(params, ns):
    # params=[nu1, nu2, t, mig]
    return moments.Demographics2D.split_mig(params, ns, pop_ids=["pop0", "pop1"])

# Returns tab-separated string encoding moments inference output
def run_moments_on_two_pop_split(sfs, sfs_params, optimization_repeats, ploidy, output_string):
    # Setup for moments
    lower_bound = [1e-3 for i in range(4)]
    upper_bound = [1 for i in range(4)]
    ns = [(int(sfs_params[0]) * ploidy) for i in range(2)]
    lls = []

    # Perform moments optimization "optimizaiton_iterations" many times in order
    # to avoid local extrema problems
    for i in range(optimization_repeats):
        params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(4)]
        popt = moments.Inference.optimize(params, sfs, two_pop_split,
                                              lower_bound=lower_bound,
                                              upper_bound=upper_bound,
                                              maxiter=10, verbose=0,
                                              multinom=True)
        model = moments.Demographics2D.split_mig(popt, ns)
        lls.append(moments.Inference.ll_multinom(model, sfs))

    # Add results to out_put string, which does not get returned because it is a global variable
    output_string += "%f\t%f" % (statistics.median(lls), max(lls))
    print(output_string)

def main(input_sfs_file, optimization_repeats, ploidy):
    input_sfs_file_name = input_sfs_file.split("/")[-1]
    # Perform moments inference on SFS created from two_pop_split demography
    if input_sfs_file_name.startswith("two_pop_split_demography"):
        # Use regex to extract numbers from string "sfs_file".
        # sfs_params = [n, msprime_seed, depth, poolseq_seed]
        sfs_params = re.findall(r'\d+', input_sfs_file_name)
        sfs = import_sfs_from_serialized_file(input_sfs_file, ["pop0", "pop1"])
        # Add metadata from name of SFS file to output_string
        output_string = "two_pop_split\t"
        for param in sfs_params:
            output_string += param + "\t"
        run_moments_on_two_pop_split(sfs, sfs_params, optimization_repeats, ploidy, output_string)

# Handle command line arguments
if len(sys.argv) < 3:
    print("Error: Two command line arguments required.")
    sys.exit()
input_sfs_file = sys.argv[1]
optimization_repeats = int(sys.argv[2])
ploidy = 2
if len(sys.argv) > 3:
    ploidy = int(sys.argv[3])

main(input_sfs_file, optimization_repeats, ploidy)

# -----
# JUNK
# -----

input_file = "/scratch/djb3ve/data/second_pipeline/sfss/two_pop_split_demography_n10_seed1_depth40_seed1.npy"
sfs = np.load(input_file)
pop_ids = ["pop0", "pop1"]
pooled_sfs_obj = moments.Spectrum(sfs, mask_corners=True, pop_ids=pop_ids)
print(pooled_sfs_obj)
print("/////")

vcf = "/scratch/djb3ve/data/second_pipeline/vcfs/two_pop_split_demography_n10_seed1.vcf"
popinfo = "/scratch/djb3ve/data/second_pipeline/popinfos/two_pop_split_demography_n10_popinfo.txt"
data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
ns = [20, 20]
initial_sfs_obj = moments.Spectrum.from_data_dict(data_dict, pop_ids=pop_ids, projections=ns, polarized=False)
print(initial_sfs_obj)
