import msprime
import sys

# Handle command line arguments to take in migration rate exponent from user input
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()
mig_rate_exp = int(sys.argv[1])

# Setup
out_dir = "/scratch/djb3ve/data/"
# Standard D. melanogaster mutation rate is 2.9e-6
mu = 1e-7
rho = 1.25e-7
seq_len = 1e5
mig_rate = 10 ** -mig_rate_exp
default_pop = 100
T_split = 5
samples = [msprime.SampleSet(10, population="pop0"),
           msprime.SampleSet(10, population="pop1")]

sample_counts = {}
for i in range(2):
    sample_counts[i] = 10

demography = msprime.Demography()
demography.add_population(name="pop0", initial_size=default_pop)
demography.add_population(name="pop1", initial_size=default_pop)
demography.add_population(name="ancestral", initial_size=default_pop * 2)
demography.add_population_split(time=T_split, derived=["pop0", "pop1"], ancestral="ancestral")
demography.set_symmetric_migration_rate(["pop0", "pop1"], mig_rate)

ts = msprime.sim_ancestry(samples=samples, demography=demography,
                          sequence_length=seq_len, recombination_rate=rho,
                          ploidy=2, random_seed=1)
mts = msprime.sim_mutations(ts, rate=mu, random_seed=1)

# Write results to VCF
with open(out_dir + "alt_2island_%dmig_model.vcf" %
          (mig_rate_exp), "w") as f:
    mts.write_vcf(f)