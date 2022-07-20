import numpy as np
import msprime
import sys
import abcpy
from abcpy.continuousmodels import Uniform
from abcpy.statistics import Identity
from abcpy.distances import Euclidean
from abcpy.backends import BackendDummy as Backend
from abcpy.inferences import PMCABC
from abcpy.probabilisticmodels import ProbabilisticModel, InputConnector
from abcpy.perturbationkernel import DefaultKernel

def generate_sfs(split_pop_size, t_split, mig_rate, sample_size=500, seq_len=1e6, rho=1.25e-12, ploidy=2, mu=2.9e-9):
    ancestral_pop_size = split_pop_size * 2
    samples={"pop0": sample_size, "pop1": sample_size}
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=ancestral_pop_size)
    dem.add_population(name="pop0", initial_size=split_pop_size)
    dem.add_population(name="pop1", initial_size=split_pop_size)
    dem.add_population_split(time=t_split, derived=["pop0", "pop1"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["pop0", "pop1"], mig_rate)
    ts = msprime.sim_ancestry(samples=samples, demography=dem,
                              sequence_length=seq_len,
                              recombination_rate=rho,
                              ploidy=ploidy)
    mts = msprime.sim_mutations(ts, rate=mu)
    return [mts.allele_frequency_spectrum(sample_sets=[mts.samples()[:sample_size * ploidy],
                                                       mts.samples()[sample_size * ploidy:]])]

class Sfs(ProbabilisticModel):
    def __init__(self, parameters, name="Sfs"):
        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        n = 10 ** input_values[0]
        t = 10 ** input_values[1]
        #m = 10 ** input_values[2]
        m = 1e-6

        return_list = [generate_sfs(n, t, m) for i in range(k)]
        return return_list

    def _check_input(self, input_values):
        return True
    def _check_output(self, values):
        return True
    def get_output_dimension(self):
        return 2

obs = generate_sfs(1e5, 1e3, 1e-6)

log_n = Uniform([[0], [7]], name="log_n")
log_t = Uniform([[0], [7]], name="log_t")
#log_m = Uniform([[-7], [0]], name="log_m")
#model = Sfs([log_n, log_t, log_m], name="model")
model = Sfs([log_n, log_t], name="model")
statistics_calculator = Identity()
distance_calculator = Euclidean(statistics_calculator)
backend = Backend()
#kernel = DefaultKernel([log_n, log_t, log_m])
kernel = DefaultKernel([log_n, log_t])
sampler = PMCABC([model], [distance_calculator], backend, kernel, seed=1)
journal = sampler.sample([obs], steps=5, epsilon_init=np.array([1]), n_samples=10)

print(journal.get_parameters())
print(journal.posterior_mean())
