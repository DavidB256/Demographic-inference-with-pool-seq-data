import msprime

out_dir = "/scratch/djb3ve/data/"
mu = 2.9e-6

dem = msprime.Demography()
dem.add_population(initial_size=100)
dem.add_instantaneous_bottleneck(time=25, strength=200, population=0)

ts = msprime.sim_ancestry(20, demography=dem, random_seed=1)
mut = ts.sim_mutations(ts, mu, random_seed=1)

# Write results to VCF
with open(out_dir + "bottleneck_model.vcf", "w") as f:
    mut.write_vcf(f)
