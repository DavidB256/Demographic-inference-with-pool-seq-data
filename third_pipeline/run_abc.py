import numpy as np
import msprime
import sys
import abcpy
from abcpy.continuousmodels import Uniform
from abcpy.statistics import Identity
from abcpy.distances import Euclidean
from abcpy.backends import BackendDummy
from abcpy.inferences import PMCABC
from abcpy.probabilisticmodels import ProbabilisticModel, InputConnector
from abcpy.perturbationkernel import DefaultKernel
import yaml
import re

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

def generate_sfs(n_1, n_2, t_split, mig_rate, sample_size,
                 seq_len=yd["dem_params"]["seq_len"],
                 rho=yd["dem_params"]["recombination_rate"],
                 ploidy=yd["dem_params"]["ploidy"],
                 mu=yd["dem_params"]["mutation_rate"]):
    ancestral_pop_size = n_1 + n_2
    samples={"pop0": sample_size, "pop1": sample_size}

    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=ancestral_pop_size)
    dem.add_population(name="pop0", initial_size=n_1)
    dem.add_population(name="pop1", initial_size=n_2)
    dem.add_population_split(time=t_split, derived=["pop0", "pop1"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["pop0", "pop1"], mig_rate)
    ts = msprime.sim_ancestry(samples=samples, demography=dem,
                              sequence_length=seq_len,
                              recombination_rate=rho,
                              ploidy=ploidy)
    mts = msprime.sim_mutations(ts, rate=mu)
    sfs = mts.allele_frequency_spectrum(sample_sets=[mts.samples()[:sample_size * ploidy],
                                                     mts.samples()[sample_size * ploidy:]])
    # Scale SFS by sequence length in order to convert its values to allele counts
    return [sfs * seq_len]

class Sfs(ProbabilisticModel):
    def __init__(self, parameters, sample_size, name="Sfs"):
        self.sample_size = sample_size
        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        n_1 = 10 ** input_values[0]
        n_2 = 10 ** input_values[1]
        t = 10 ** input_values[2]
        m = 10 ** input_values[3]

        return_list = [generate_sfs(n_1, n_2, t, m, self.sample_size) for i in range(k)]
        return return_list

    def _check_input(self, input_values):
        return True
    def _check_output(self, values):
        return True
    def get_output_dimension(self):
        return 2

if __name__ == "__main__":
    # Handle command line arguments
    sfs_name = sys.argv[1]
    # Use regex to extract numbers "sfs_name" into a list
    sfs_name_numbers = re.findall(r'\d+', sfs_name)
    observed_sfs_file = yd["pipeline_params"]["data_dir"] + "sfss/" + sfs_name
    sample_size = int(sfs_name_numbers[0])
    obs = [np.load(observed_sfs_file)]

    # Parameters to be estimated by ABC
    log_n_1 = Uniform(yd["abc_params"]["prior_bounds"]["log_n_1"], name="log_n_1")
    log_n_2 = Uniform(yd["abc_params"]["prior_bounds"]["log_n_2"], name="log_n_2")
    log_t = Uniform(yd["abc_params"]["prior_bounds"]["log_t_split"], name="log_t")
    log_m = Uniform(yd["abc_params"]["prior_bounds"]["log_m"], name="log_m")
    param_list = [log_n_1, log_n_2, log_t, log_m]

    # Set up for ABC run
    model = Sfs(param_list, sample_size, name="model")
    statistics_calculator = Identity()
    distance_calculator = Euclidean(statistics_calculator)
    backend = BackendDummy()
    kernel = DefaultKernel(param_list)
    sampler = PMCABC([model], [distance_calculator], backend, kernel, seed=1)

    print("Starting ABC run.")

    # Run ABC
    journal = sampler.sample([obs], steps=5, epsilon_init=np.array([1e9]), n_samples=100)

    # Extract statistics from journal object
    means = journal.posterior_mean()
    covariance = journal.posterior_cov()[0]

    # Write to output file
    output_list = sfs_name_numbers
    for mean in means:
        output_list.append(str(means[mean]))
    for i in range(len(covariance)):
        output_list.append(str(covariance[i][i]))
    with open(yd["pipeline_params"]["data_dir"] + yd["pipeline_params"]["output_file_name"], "a") as f:
        for output in output_list:
            f.write(output + "\t")
        f.write("\n")
















