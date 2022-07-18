import numpy as np
import msprime
import abcpy
from abcpy.continuousmodels import Uniform
from abcpy.statistics import Identity
from abcpy.distances import Euclidean, Wasserstein
from abcpy.perturbationkernel import DefaultKernel
from abcpy.backends import BackendDummy as Backend
from abcpy.inferences import RejectionABC
from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector

def generate_sfs(n, samples={"pop0": 10}, seq_len=5e5, rho=1.25e-12, ploidy=2, mu=2.9e-9):
    dem = msprime.Demography()
    dem.add_population(name="pop0", initial_size=n)
    ts = msprime.sim_ancestry(samples=samples, demography=dem,
                              sequence_length=seq_len,
                              recombination_rate=rho,
                              ploidy=ploidy)
    mts = msprime.sim_mutations(ts, rate=mu)
    return [mts.allele_frequency_spectrum()]

class Sfs(ProbabilisticModel, Continuous):
    def __init__(self, parameters, name="Sfs"):
        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        n = 10 ** input_values[0]

        return_list = [generate_sfs(n) for i in range(k)]
        return return_list

    def _check_input(self, input_values):
        return True
    def _check_output(self, values):
        return True
    def get_output_dimension(self):
        return 2

obs = generate_sfs(1000)

log_n = Uniform([[1e-3], [5]], name="log_n")
model = Sfs([log_n], name="model")
statistics_calculator = Identity()
distance_calculator = Euclidean(statistics_calculator)
backend = Backend()
sampler = RejectionABC([model], [distance_calculator], backend, seed=1)
n_samples = 5
n_samples_per_param = 1
journal = sampler.sample([obs], n_samples=n_samples, n_samples_per_param=n_samples_per_param,
                         epsilon=2e-5)

print(journal.get_parameters())
avg = 0
for est in journal.get_parameters()["log_n"]:
    avg += est[0][0]
print(avg / n_samples)
