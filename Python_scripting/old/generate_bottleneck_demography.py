import msprime

out_dir = "/scratch/djb3ve/data/"
mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e5

dem = msprime.Demography()
dem.add_population(initial_size=100)
dem.add_instantaneous_bottleneck(time=25, strength=200, population=0)

ts = msprime.sim_ancestry(samples=5, demography=dem,
                          sequence_length=seq_len, random_seed=1,
                          recombination_rate=rho)
mut = msprime.sim_mutations(ts, mu, random_seed=1)

# Write results to VCF
with open(out_dir + "bottleneck_model.vcf", "w") as f:
    mut.write_vcf(f)
