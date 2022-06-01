import msprime

# Setup
out_dir = "/scratch/djb3ve/data/"
mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e5
mig_rate = 0.001
default_pop = 100

sample_counts = {}
for i in range(3):
    sample_counts[i] = 2

# Prepare demography
dem = msprime.Demography.island_model([default_pop for i in range(3)],
                                      migration_rate=mig_rate)

# Perform ancestry simulation
ts = msprime.sim_ancestry(samples=sample_counts, demography=dem,
                          sequence_length=seq_len, random_seed=1,
                          recombination_rate=rho)
# Add mutations
mut = msprime.sim_mutations(ts, rate=mu, random_seed=1)

# Write results to VCF
with open(out_dir + "3island_small_model.vcf", "w") as f:
    mut.write_vcf(f)
