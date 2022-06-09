import msprime
import sys
import os.path

# Setup
overwrite = False
if len(sys.argv) > 1 and sys.argv[1] in ["o", "overwrite"]:
    overwrite = True

output_dir = "/scratch/djb3ve/data/first_models/"

mu = 2.9e-6 # D. melanogaster mutation rate
rho = 1.25e-7 # D. melanogaster recombination rate
seq_len = 1e5
pop_size = 100
mig_rate = 1e-2
t_split = 10

sample_sizes = [10, 50, 100, 200, 500]

def write_vcf_from_demography(file_name, samples, demography,
                              sequence_length=seq_len, recombination_rate=rho,
                              ploidy=2, random_seed=1, mutation_rate=mu):
    ts = msprime.sim_ancestry(samples=samples, demography=demography,
                              sequence_length=sequence_length,
                              recombination_rate=recombination_rate,
                              ploidy=ploidy, random_seed=random_seed)
    mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=random_seed)
    with open(output_dir + file_name, "w") as f:
        mts.write_vcf(f)


for sample_size in sample_sizes:
    # control
    file_name = f"control_demography_{sample_size}n.vcf"
    if not os.path.exists(output_dir + file_name) or overwrite:
        demography = msprime.Demography()
        demography.add_population(name="pop0", initial_size=pop_size)
        write_vcf_from_demography(file_name, sample_size, demography)

    # two_pop_split
    file_name = f"two_pop_split_{sample_size}n.vcf"
    if not os.path.exists(output_dir + file_name) or overwrite:
        demography = msprime.Demography()
        demography.add_population(name="pop0", initial_size=pop_size)
        demography.add_population(name="pop1", initial_size=pop_size)
        demography.add_population(name="ancestral", initial_size=pop_size * 2)
        demography.add_population_split(time=t_split, derived=["pop0", "pop1"], ancestral="ancestral")
        demography.set_symmetric_migration_rate(["pop0", "pop1"], mig_rate)
        write_vcf_from_demography(file_name,
                                  {"pop0": sample_size, "pop1": sample_size},
                                  demography)

    # bottleneck
    # file_name = f"two_pop_split_{sample_size}n.vcf"
    # if not os.path.exists(output_dir + file_name) or overwrite:
    #     demography = msprime.Demography()
    #     demography.add_population
