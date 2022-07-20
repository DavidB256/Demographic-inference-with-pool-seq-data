import numpy as np
import msprime
import sys
import abcpy
from abcpy.continuousmodels import Uniform
from abcpy.statistics import Identity
from abcpy.distances import Euclidean
from abcpy.backends import BackendMPI as Backend
from abcpy.inferences import RejectionABC
from abcpy.probabilisticmodels import ProbabilisticModel, InputConnector

def generate_sfs(split_pop_size, t_split=1000, sample_size=10,
                 seq_len=1e6, rho=1.25e-12, ploidy=2, mu=2.9e-9, mig_rate=1e-5):
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
    sfs = mts.allele_frequency_spectrum(sample_sets=[mts.samples()[:sample_size * ploidy],
                                                     mts.samples()[sample_size * ploidy:]])
    return sfs[1][1]

class Sfs(ProbabilisticModel):

    def __init__(self, parameters, name='Sfs'):
        # We expect input of type parameters = [mu, sigma]
        if not isinstance(parameters, list):
            raise TypeError('Input of Normal model is of type list')

        if len(parameters) != 1:
            raise RuntimeError('Input list must be of length 1, containing [log_n].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 1:
            raise ValueError('Number of parameters are 1.')
        return True

    def _check_output(self, values):
        if not isinstance(values, np.ndarray):
            raise ValueError('Output of the normal distribution is always a numpy array.')
        if values.shape[0] != 2:
            raise ValueError('Output shape should be of dimension 2.')
        return True

    def get_output_dimension(self):
        return 2

    def forward_simulate(self, input_values, k, rng=np.random.RandomState, mpi_comm=None):
        if mpi_comm is None:
            ValueError('MPI-parallelized simulator model needs to have access \
            to a MPI communicator object')
        print("Start Forward Simulate on rank {}".format(mpi_comm.Get_rank()))
        rank = mpi_comm.Get_rank()
        # Extract the input parameters
        log_n = input_values[rank]
        n = 10 ** log_n
        # Do the actual forward simulation
        vector_of_k_samples = np.array([generate_sfs(n) for i in range(k)])
        print(f"vector_of_k_samples: {vector_of_k_samples}")

        # Send everything back to rank 0
        data = mpi_comm.gather(vector_of_k_samples, root=0)
        print(f"data: {data}")

        # Format the output to obey API and broadcast it before return
        result = None
        if rank == 0:
            result = [None] * k
            for i in range(k):
                element0 = data[0][i]
                element1 = data[1][i]
                point = np.array([element0, element1])
                result[i] = point
            result = [np.array([result[i]]).reshape(-1, ) for i in range(k)]
            print(f"End forward sim on master {result}")
            return result
        else:
            print("End forward sim on workers")
            return None

obs = [np.array([generate_sfs(10000)])]

log_n = Uniform([[0], [5]], name="log_n")
model = Sfs([log_n], name="model")
statistics_calculator = Identity(degree=2, cross=False)
distance_calculator = Euclidean(statistics_calculator)
backend = Backend(process_per_model=2)
sampler = RejectionABC([model], [distance_calculator], backend, seed=1)
epsilon = float(sys.argv[1])
print(obs)
journal = sampler.sample([obs], n_samples=2, n_samples_per_param=1, epsilon=1)

print(journal.get_parameters())
print(journal.posterior_mean())
