import msprime

# Generates small VCF file representing a single isolated subpopulation for use
# testing purposes

# Setup
mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e5
default_pop = 100
out_dir = "/scratch/djb3ve/data/"

# Perform simulation
dem = msprime.Demography.isolated_model([default_pop])
ts = msprime.sim_ancestry(samples=5, demography=dem,
                          sequence_length=seq_len, random_seed=1,
                          recombination_rate=rho)
mut = msprime.sim_mutations(ts, rate=mu, random_seed=1)

# Write to VCF
with open(out_dir + "small_vcf.vcf", "w") as f:
    mut.write_vcf(f)
