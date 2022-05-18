import msprime

mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e6
sample_counts = {0: 100, 1: 100, 2: 100, 3: 100}
mig_rate = 0.1
default_pop = 100
out_dir = "/scratch/djb3ve/msprime/may16/dir_model_vcfs/"

dem_control = msprime.Demography.isolated_model([default_pop])
dem_continent_islands = msprime.Demography.island_model([default_pop * 100, default_pop, default_pop, default_pop], mig_rate)
dem_equal_islands = msprime.Demography.island_model([default_pop] * 4, mig_rate)
dem_stepping_stone = msprime.Demography.stepping_stone_model([default_pop] * 4, mig_rate)

ts_control = msprime.sim_ancestry(samples=100, demography=dem_control, sequence_length=seq_len, random_seed=1, recombination_rate=rho)
ts_continent_islands = msprime.sim_ancestry(samples=sample_counts, demography=dem_continent_islands, sequence_length=seq_len, random_seed=1, recombination_rate=rho)
ts_equal_islands = msprime.sim_ancestry(samples=sample_counts, demography=dem_equal_islands, sequence_length=seq_len, random_seed=1, recombination_rate=rho)
ts_stepping_stone = msprime.sim_ancestry(samples=sample_counts, demography=dem_stepping_stone, sequence_length=seq_len, random_seed=1, recombination_rate=rho)

mut_control = msprime.sim_mutations(ts_control, rate=mu, random_seed=1)
mut_continent_islands = msprime.sim_mutations(ts_continent_islands, rate=mu, random_seed=1)
mut_equal_islands = msprime.sim_mutations(ts_equal_islands, rate=mu, random_seed=1)
mut_stepping_stone = msprime.sim_mutations(ts_stepping_stone, rate=mu, random_seed=1)

with open(out_dir + "control.vcf", "w") as f:
    mut_control.write_vcf(f)
with open(out_dir + "continent_islands.vcf", "w") as f:
    mut_continent_islands.write_vcf(f)
with open(out_dir + "equal_islands.vcf", "w") as f:
    mut_equal_islands.write_vcf(f)
with open(out_dir + "stepping_stone.vcf", "w") as f:
    mut_stepping_stone.write_vcf(f)
