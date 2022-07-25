import sys
import os
import numpy as np
import pandas as pd
import yaml
import re
import msprime

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

# Returns number of non-reference alleles in list argument "row", whose elements are
# pipe-separated VCF allele records.
# Treats all sites as biallelic, i.e. all alternate alleles are treated as a single
# "non-reference" allele.
def get_allele_count_from_vcf_row(row):
    cumulative_count = 0
    for polymorphism in row:
        allele_list = polymorphism.split("|")
        alt_allele_count = len(allele_list) - allele_list.count("0")
        cumulative_count += alt_allele_count
    return cumulative_count

# Replicate of "sample.alleles(mode="coverage")" in Thomas Taus' poolSeq
# (https://github.com/ThomasTaus/poolSeq/blob/master/R/simaf.R).
# Lacks support for specifying pool-seq coverage at each site with vector argument for "size" parameter.
# Lacks support for mode="individuals".
def get_pooled_allele_freqs(allele_freqs, poolseq_depth):
    cov = np.random.poisson(lam=poolseq_depth, size=len(allele_freqs))
    cov = np.where(cov == 0, 1, cov)
    p_sampled = np.random.binomial(n=cov, p=allele_freqs, size=len(allele_freqs)) / cov
    return p_sampled

# Replicate of moments' "Spectrum.from_data_dict" function, but with the
# application of pool-seq noise to allele frequencies
# (https://moments.readthedocs.io/en/devel/api/api_moments.html#moments.Spectrum.from_data_dict).
# The rounding methods are taken from genomalicious' "dadi_inputs_pools" function
# (https://rdrr.io/github/j-a-thia/genomalicious/man/dadi_inputs_pools.html).
def get_folded_sfs(vcf_file, popinfo_file, poolseq_depth, haploid_counts):
    # Load in files as pandas data frames
    vcf_df = pd.read_table(vcf_file, header=5)
    # Remove multiallelic polymorphisms, i.e. those with more than one alternate allele
    if yd["pipeline_params"]["ignore_multiallelic_site"]:
        vcf_df = vcf_df[vcf_df.ALT.apply(lambda x: len(x) == 1)]
        vcf_df.reset_index(drop=True, inplace=True)
    # Load popinfo file into dataframe
    popinfo_df = pd.read_table(popinfo_file, header=None)

    # Number of rows of "vcf_df" is the number of polymorphic sites recorded in the VCF
    num_of_polymorphisms = vcf_df.shape[0]

    # Subdivide columns of VCF file by population, with one population per element of
    # "vcf_dfs_by_population"
    vcf_dfs_by_population = []
    for pop_id in popinfo_df[1].unique():
        tsks = popinfo_df[popinfo_df[1] == pop_id][0]
        vcf_dfs_by_population.append(vcf_df[tsks])

    allele_counts_by_population = []
    for pop_vcf_df in vcf_dfs_by_population:
        allele_counts = pop_vcf_df.apply(get_allele_count_from_vcf_row, axis=1)
        allele_counts_by_population.append(allele_counts)

    # If "poolseq_depth" is nonzero, then applying pool-seq noise.
    if poolseq_depth:
        # Convert population-based subsets of VCF rows into allele frequencies
        allele_freqs_by_population = []
        for i, pop_allele_counts in enumerate(allele_counts_by_population):
            allele_freqs_by_population.append(pop_allele_counts / haploid_counts[i])

        # Apply pool-seq noise to allele frequencies
        poolseq_allele_freqs_by_population = []
        for pop_allele_freqs in allele_freqs_by_population:
            poolseq_allele_freqs = get_pooled_allele_freqs(pop_allele_freqs, poolseq_depth)
            poolseq_allele_freqs_by_population.append(poolseq_allele_freqs)

        # Use a rounding method adapted from genomalicious to obtain allele counts from frequencies
        allele_counts_by_population = pd.DataFrame()
        # The "counts" method intuitively multiplies pooled allele frequencies by haploid counts
        # and then rounds to the nearest integer to yield allele counts. In order to avoid
        # erasing rare alleles, pre-rounding allele counts less than one are always rounded up.
        if yd["pipeline_params"]["rounding_method"] == "counts":
            for i, pop_pooled_allele_freqs in enumerate(poolseq_allele_freqs_by_population):
                # Scale frequences by haploid counts to achieve allele counts once rounding has been performed
                pop_pooled_allele_counts = pop_pooled_allele_freqs * haploid_counts[i]
                # Prevent rare alleles from getting rounded down to zero, and thus removed from the model
                for j, count in enumerate(pop_pooled_allele_counts):
                    if 0 < count < 0.5:
                        pop_pooled_allele_counts[j] = 1
                # Round scaled allele frequencies to achieve allele counts
                pop_pooled_allele_counts = np.round(pop_pooled_allele_counts)
                allele_counts_by_population[i] = pop_pooled_allele_counts
        # The "probs" method uses binomial sampling in order to avoid the issue of having to round
        # rare alleles counts up in the "counts" method. This introduces more noise and
        # is promising, but should probably not be used as a default.
        # THIS HAS NOT BEEN TESTED YET.
        elif yd["pipeline_params"]["rounding_method"] == "probs":
            for i, pop_pooled_allele_freqs_and_coverage in enumerate(poolseq_allele_freqs_by_population):
                allele_counts_by_population[i] = np.random.binomial(n=haploid_counts[i],
                                                                    p=pop_pooled_allele_freqs,
                                                                    size=num_of_polymorphisms)
        else:
            print("Error: Invalid rounding method.")
            sys.exit()
    else:
        allele_counts_by_population = pd.DataFrame(allele_counts_by_population).transpose()

    # Fold allele counts so that the minor allele is counted as the alternate allele at each site.
    # I am unsure of how SFSs are constructed for polyploidal populations, so this might not work for
    # non-diploid simulations.
    total_haploids_sampled = np.sum(haploid_counts)
    def fold_allele_counts(row):
        if np.sum(row) >= total_haploids_sampled / 2:
            for i, count in enumerate(row):
                #row[i] = min(count, haploid_counts[i] - count)
                row[i] = haploid_counts[i] - count
        return row
    folded_allele_counts_by_population = allele_counts_by_population.apply(fold_allele_counts, axis=1)

    # Construct numpy array containing SFS from pooled_allele_counts_by_population
    incremented_haploid_counts = [i + 1 for i in haploid_counts]
    sfs = np.zeros(incremented_haploid_counts)
    for i, row in folded_allele_counts_by_population.iterrows():
        sfs[tuple(row.astype(int))] += 1

    # Mask corners, which represent alleles that are absent and ubiquitous, respectively
    absent_allele_coords = [0 for i in haploid_counts]
    sfs[tuple(absent_allele_coords)] = 0
    sfs[tuple(haploid_counts)] = 0

    return sfs

def generate_ts(sample_size):
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
    return mts

def print_sfss(sfs1, sfs2):
    sfs = sfs1 - sfs2
    for i in sfs:
        for j in i:
            print(j, end=" ")
        print()
    print()

if __name__ == "__main__":
    # Prevent numpy from using scientific notation
    np.set_printoptions(suppress=True)

    vcf_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/qa_testing/vcf.vcf"
    output_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/qa_testing/sfs_distances_seqlen6.txt"

    #sample_sizes = [10, 50, 100, 200] * 40
    sample_sizes = [10] * 40

    with open(output_file, "w") as g:
        for sample_size in sample_sizes:
            print(sample_size)

            ts = generate_ts(sample_size)
            with open(vcf_file, "w") as f:
                ts.write_vcf(f)

            haploids = sample_size * yd["dem_params"]["ploidy"]
            msprime_sfs = ts.allele_frequency_spectrum(sample_sets=[ts.samples()[:haploids], ts.samples()[haploids:]])
            msprime_sfs *= yd["dem_params"]["seq_len"]

            popinfo_file = f"/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data_v1/popinfos/n{sample_size}_popinfo.txt"
            haploid_counts = [haploids] * 2
            my_sfs = get_folded_sfs(vcf_file, popinfo_file, 0, haploid_counts)

            dist = np.linalg.norm(msprime_sfs - my_sfs)

            if dist > sample_size:
                print_sfss(msprime_sfs, my_sfs)

            g.write(str(sample_size) + "\t" + str(dist) + "\n")



















