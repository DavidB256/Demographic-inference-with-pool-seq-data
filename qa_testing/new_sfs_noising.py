import numpy as np
import yaml
import msprime

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

def print_sfs(sfs):
    for i in sfs:
        for j in i:
            print(j, end=" ")
        print()
    print()

def add_noise_to_sfs(sfs, poolseq_depth, method="counts"):
    haploid_counts = [i + 1 for i in sfs.shape]
    noised_sfs = np.zeros(sfs.shape)
    it  = np.nditer(sfs, flags=['multi_index'])

    for i in it:
        freqs = [index / haploid_counts[j] for j, index in enumerate(it.multi_index)]
        for allele in range(i):
            # Generate an allele's frequencies in each subpopulation
            coverages = np.random.poisson(lam=poolseq_depth, size=sfs.ndim)
            coverages = np.where(coverages == 0, 1, coverages)
            pooled_freqs = np.random.binomial(n=coverages, p=freqs, size=sfs.ndim) / coverages

            # Convert frequencies to counts via a rounding method
            if method == "counts":
                pooled_counts = np.multiply(pooled_freqs, haploid_counts)
                for j, freq in enumerate(pooled_counts):
                    if 0 < freq < 0.5:
                        pooled_counts[j] = 1
                    pooled_counts = np.round(pooled_counts).astype(int)
            elif method == "probs":
                pooled_counts = np.random.binomial(n=haploid_counts, p=pooled_freqs, size=sfs.ndim)
            else:
                raise "Unrecognized rounding method."

            # Fold counts
            if np.sum(pooled_counts) > np.sum(haploid_counts) / 2:
                for j, count in enumerate(pooled_counts):
                    pooled_counts[j] = min(pooled_counts[j],
                                           haploid_counts[j] - pooled_counts[j])

            # Mask corners
            absent_allele_coords = [0] * sfs.ndim
            noised_sfs[tuple(absent_allele_coords)] = 0
            noised_sfs[tuple(haploid_counts)] = 0

            # Increment corresponding element
            noised_sfs[tuple(pooled_counts)] += 1

    return noised_sfs.astype(int)

def generate_sfs_from_msprime(sample_size):
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=yd["dem_params"]["ancestral_population_size"])
    dem.add_population(name="pop0", initial_size=yd["dem_params"]["split_population_size"])
    dem.add_population(name="pop1", initial_size=yd["dem_params"]["split_population_size"])
    dem.add_population_split(time=yd["dem_params"]["t_split"], derived=["pop0", "pop1"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["pop0", "pop1"], yd["dem_params"]["mig_rate"])

    samples = {"pop0": sample_size, "pop1": sample_size}
    ts = msprime.sim_ancestry(samples=samples, demography=dem,
                              sequence_length=yd["dem_params"]["seq_len"],
                              recombination_rate=yd["dem_params"]["recombination_rate"],
                              ploidy=yd["dem_params"]["ploidy"])
    mts = msprime.sim_mutations(ts, rate=yd["dem_params"]["mutation_rate"])

    haploids = sample_size * yd["dem_params"]["ploidy"]
    sfs = mts.allele_frequency_spectrum(sample_sets=[ts.samples()[:haploids],
                                                     ts.samples()[haploids:]])
    return (sfs * yd["dem_params"]["seq_len"]).astype(int)

sample_size = 10
output_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/qa_testing/new_sfs_noising_dists_probs.txt"

with open(output_file, "a") as f:
    for depth in yd["pipeline_params"]["poolseq_depths"] + [1000]:
        print("Starting run with depth %d." % depth)
        for i in range(40):
            print("Starting iteration %d." % i)
            sfs = generate_sfs_from_msprime(sample_size)
            noised_sfs = add_noise_to_sfs(sfs, depth, "probs")

            dist = np.linalg.norm(sfs - noised_sfs)

            f.write("%d\t%f\n" % (depth, dist))
            print("%d\t%f\n" % (depth, dist))












