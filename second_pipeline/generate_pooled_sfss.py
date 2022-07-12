import sys
import os
import numpy as np
import pandas as pd

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
def get_pooled_allele_freqs(allele_freqs, poolseq_depth, ploidy=2):
    cov = np.random.poisson(lam=poolseq_depth, size=len(allele_freqs))
    cov = np.where(cov == 0, 1, cov)
    p_sampled = np.random.binomial(n=cov, p=allele_freqs, size=len(allele_freqs)) / cov
    return p_sampled

# Replicate of moments' "Spectrum.from_data_dict" function, but with the
# application of pool-seq noise to allele frequencies
# (https://moments.readthedocs.io/en/devel/api/api_moments.html#moments.Spectrum.from_data_dict).
# The rounding methods are taken from genomalicious' "dadi_inputs_pools" function
# (https://rdrr.io/github/j-a-thia/genomalicious/man/dadi_inputs_pools.html).
def get_pooled_folded_fs(vcf_file, popinfo_file, haploid_counts, poolseq_depth, rounding_method="counts", ploidy=2):
    # Load in files as pandas data frames
    vcf_df = pd.read_table(vcf_file, header=5)
    # Remove non-biallelic polymorphisms
    vcf_df = vcf_df.loc[vcf_df["ALT"].isin(["A", "T", "G", "C"])]
    popinfo_df = pd.read_table(popinfo_file, header=None)

    # Number of rows of "vcf_df" is the number of polymorphic sites recorded in the VCF
    num_of_polymorphisms = vcf_df.shape[0]

    # Subdivide columns of VCF file by population, with one population per element of
    # "vcf_dfs_by_population"
    vcf_dfs_by_population = []
    for pop_id in popinfo_df[1].unique():
        tsks = popinfo_df[popinfo_df[1] == pop_id][0]
        vcf_dfs_by_population.append(vcf_df[tsks])

    allele_freqs_by_population = []
    for i, pop_vcf_df in enumerate(vcf_dfs_by_population):
        allele_counts = pop_vcf_df.apply(get_allele_count_from_vcf_row, axis=1)
        allele_freqs_by_population.append(allele_counts / haploid_counts[i])

    poolseq_allele_freqs_and_coverage_by_population = []
    for pop_allele_freqs in allele_freqs_by_population:
        poolseq_allele_freqs_and_coverage = get_pooled_allele_freqs(pop_allele_freqs, poolseq_depth, ploidy)
        poolseq_allele_freqs_and_coverage_by_population.append(poolseq_allele_freqs_and_coverage)

    pooled_allele_counts_by_population = []
    # The "counts" method intuitively multiplies pooled allele frequencies by haploid counts
    # and then rounds to the nearest integer to yield allele counts. In order to avoid
    # erasing rare alleles, pre-rounding allele counts less than one are always rounded up.
    if rounding_method == "counts":
        for i, pop_pooled_allele_freqs_and_coverage in enumerate(poolseq_allele_freqs_and_coverage_by_population):
            # Scale frequences by haploid counts to achieve allele counts once rounding has been performed
            pop_pooled_allele_counts = pop_pooled_allele_freqs_and_coverage * haploid_counts[i]
            # Prevent rare alleles from getting rounded down to zero, and thus removed from the model
            for j, count in enumerate(pop_pooled_allele_counts):
                if 0 < count < 0.5:
                    pop_pooled_allele_counts[j] = 1
            # Round scaled allele frequencies to achieve allele counts
            pop_pooled_allele_counts = np.round(pop_pooled_allele_counts)
            pooled_allele_counts_by_population.append(pop_pooled_allele_counts)
    # The "probs" method uses binomial sampling in order to avoid the issue of having to round
    # rare alleles counts up in the "counts" method. This introduces more noise and
    # is promising, but should probably not be used as a default.
    # THIS HAS NOT BEEN TESTED YET.
    elif rounding_method == "probs":
        for i, pop_pooled_allele_freqs_and_coverage in enumerate(poolseq_allele_freqs_and_coverage_by_population):
            pop_pooled_allele_counts = np.random.binomial(n=haploid_counts[i],
                                                          p=pop_pooled_allele_freqs_and_coverage,
                                                          size=num_of_polymorphisms)
            pooled_allele_counts_by_population.append(pop_pooled_allele_counts)
    else:
        print("Error: Invalid rounding method.")
        sys.exit()

    # Construct numpy array containing SFS from pooled_allele_counts_by_population
    incremented_haploid_counts = [i + 1 for i in haploid_counts]
    sfs = np.zeros(incremented_haploid_counts)
    for i in range(num_of_polymorphisms):
        coords = []
        for pop_pooled_allele_counts in pooled_allele_counts_by_population:
            coords.append(int(pop_pooled_allele_counts[i]))
        sfs[tuple(coords)] += 1

    return sfs

def main(vcf_file, popinfo_file, haploid_counts, poolseq_depth, output_file_without_suffix, num_of_replicates, rounding_method="counts", ploidy=2):
    # Iterate once per replicate, incrementing the seed between replicates
    for seed in range(1, num_of_replicates+1):
        np.random.seed(seed)
        # Get and serialize SFS
        sfs = get_pooled_folded_fs(vcf_file, popinfo_file, haploid_counts, poolseq_depth, rounding_method, ploidy)
        np.save(f"{output_file_without_suffix}_seed{seed}.npy", sfs)

# Prevent numpy from using scientific notation
np.set_printoptions(suppress=True)

# Set up working environment
working_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"
for subdir in ["vcfs/", "serialized_pooled_sfss/"]:
    if not os.path.exists(working_dir + subdir):
        os.makedirs(working_dir + subdir)

# Handle command line arguments, which will be input from "pipeline_instructions.txt",
# except for "num_of_replicates", which is added to the "$OPTS" variable in the wrapper
# script.
vcf_file = working_dir + "vcfs/" + sys.argv[1]
popinfo_file = working_dir + "popinfos/" + sys.argv[2]
output_file_without_suffix = working_dir + "sfss/" + sys.argv[1][:-4] + "_depth" + sys.argv[4]
haploid_counts = eval(sys.argv[3])
poolseq_depth = int(sys.argv[4])
rounding_method = sys.argv[5]
ploidy = int(sys.argv[6])
num_of_replicates = int(sys.argv[7])


main(vcf_file, popinfo_file, haploid_counts, poolseq_depth, output_file_without_suffix, num_of_replicates, rounding_method, ploidy)
