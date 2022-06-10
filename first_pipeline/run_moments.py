import moments
import numpy as np

def control_model(params, ns):
    # params = [nu, T]
    return moments.Demographics1D.growth(params, ns)

def two_pop_split_model():
    pass

# This script should iterate through pipeline_instructions
instructions = "/scratch/djb3ve/data/pipeline_instructions.txt"

vcf = "/scratch/djb3ve/data/first_models/control_demography_n10_seed1.vcf"
popinfo = "/scratch/djb3ve/first_models/control_demography_n10_popinfo.txt"
ns = [10]

data_dict = moments.Misc.make_data_dict_vcf(vcf, popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"], projections=ns, polarized=False)

lower_bound = [1e-3 for i in range(2)]
upper_bound = [10 for i in range(2)]

params = [np.random.uniform(lower_bound[j], upper_bound[j]) for j in range(2)]
popt = moments.Inference.optimize_log(params, fs, control_model,
                                      lower_bound=lower_bound,
                                      upper_bound=upper_bound,
                                      maxiter=100, verbose=1)
model = control_model(popt, ns)
ll_model = moments.Inference.ll_multinom(model, fs)
theta = moments.Inference.optimal_sfs_scaling(model, fs)
