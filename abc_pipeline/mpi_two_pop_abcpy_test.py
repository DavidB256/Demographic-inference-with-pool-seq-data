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

def generate_sfs(t_split, samples={"pop0": 10, "pop1": 10}, seq_len=5e5, rho=1.25e-12, ploidy=2, mu=2.9e-9,
                 ancestral_pop_size=2e5, split_pop_size=1e5, mig_rate=1e-5):
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
    return  mts.allele_frequency_spectrum()[1]

class Sfs(ProbabilisticModel):
    def __init__(self, parameters, name="Sfs"):
        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def forward_simulate(self, input_values, k, rng=np.random.RandomState(), mpi_comm=None):
        if mpi_comm is None:
            ValueError('MPI-parallelized simulator model needs to have access \
            to a MPI communicator object')
        print("Start Forward Simulate on rank {}".format(mpi_comm.Get_rank()))
        rank = mpi_comm.Get_rank()
        # Extract the input parameters
        n = 10 ** input_values[rank]

        # Perform forward simulation
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
            # print("End forward sim on workers")
            return None

    def _check_input(self, input_values):
        return True
    def _check_output(self, values):
        return True
    def get_output_dimension(self):
        return 2

obs = [np.array([generate_sfs(1000)])]

log_t = Uniform([[0], [5]], name="log_t")
model = Sfs([log_t], name="model")
statistics_calculator = Identity(degree=2, cross=False)
distance_calculator = Euclidean(statistics_calculator)
backend = Backend(process_per_model=2)
sampler = RejectionABC([model], [distance_calculator], backend, seed=1)
epsilon = float(sys.argv[1])
journal = sampler.sample([obs], n_samples=2, n_samples_per_param=1, epsilon=epsilon)

print(journal.get_parameters())
print(journal.posterior_mean())
