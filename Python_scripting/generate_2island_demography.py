import msprime
import sys

out_dir = "/scratch/djb3ve/data/"

# Handle command line arguments
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()
mig_rate_exp = int(sys.argv[1])

# Setup
mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e5
mig_rate = 10 ** -mig_rate_exp
default_pop = 100

sample_counts = {}
for i in range(2):
    sample_counts[i] = 10

# Prepare demography
dem = msprime.Demography.island_model([default_pop for i in range(2)],
                                      migration_rate=mig_rate)

# Perform simulations
ts = msprime.sim_ancestry(samples=sample_counts, demography=dem,
                          sequence_length=seq_len, random_seed=1,
                          recombination_rate=rho)
mut = msprime.sim_mutations(ts, rate=mu, random_seed=1)

# Write results to VCF
with open(out_dir + "2island_%dmig_model.vcf" %
          (mig_rate_exp), "w") as f:
    mut.write_vcf(f)
