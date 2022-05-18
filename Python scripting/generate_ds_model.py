import msprime
import sys

# Requires one command line argument, a positive integer for s, the side length of the grid of subpopulations
# Uses msprime's island_model demography to replicate the D_s model in “Interpreting principal component analyses of spatial population genetic variation” (Novembre and Stephens, 2008)

# Handle command line arguments
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()

# Setup
mu = 2.9e-6
rho = 1.25e-7
seq_len = 1e6
mig_rate = 0.1
default_pop = 100
out_dir = "/scratch/djb3ve/msprime/may17/"
dim = sys.argv[1]
dim2 = dim ** 2

sample_counts = {}
for i in range(dim2):
    sample_counts[i] = 10

# Prepare demography
dem = msprime.Demography.island_model([default_pop for i in range(dim2)], migration_rate=0)
for i in range(dim2):
    for j in range(dim2):
        if [abs((i//dim) - (j//dim)), abs((i % dim) - (j % dim))] in [[0, 1], [1, 0]]:
            dem.set_migration_rate(i, j, mig_rate)

# Perform simulations
ts = msprime.sim_ancestry(samples=sample_counts, demography=dem, sequence_length=seq_len, random_seed=1, recombination_rate=rho)
mut = msprime.sim_mutations(ts, rate=mu, random_seed=1)

# Write results to VCF
with open("d%d_model.vcf" % (dim), "w") as f:
    mut.write_vcf(f)
