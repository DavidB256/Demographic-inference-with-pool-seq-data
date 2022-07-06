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
        alt_allele_count = len(allele_list) - allele_list.count(0)
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
    return np.column_stack(p_sampled, cov)

# Replicate of moments' "Spectrum.from_data_dict" function, but with the
# application of pool-seq noise to allele frequencies
# (https://moments.readthedocs.io/en/devel/api/api_moments.html#moments.Spectrum.from_data_dict).
# The rounding methods are taken from genomalicious' "dadi_inputs_pools" function
# (https://rdrr.io/github/j-a-thia/genomalicious/man/dadi_inputs_pools.html).
def get_pooled_folded_fs(vcf_file, popinfo_file, haploid_counts, poolseq_depth, rounding_method="counts", ploidy=2):
    # Load in files as pandas data frames
    vcf_df = pd.read_table(vcf_file, header=5)
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

    pooled_allele_freqs_by_population = []
    for pop_allele_freqs in allele_freqs_by_population:
        pooled_allele_freqs = get_pooled_allele_freqs(pop_allele_freqs, poolseq_depth, ploidy)
        pooled_allele_freqs_by_population.append(pooled_allele_freqs)

    pooled_allele_counts_by_population = []
    for i, pop_pooled_allele_freqs in enumerate(pooled_allele_freqs_by_population):
        if rounding_method == "counts":
            pooled_allele_counts = pop_pooled_allele_freqs * haploid_counts[i]
            pooled_allele_counts = np.where(pooled_allele_counts > 0 && pooled_allele_counts < 1,
                                            1, pooled_allele_counts)
        elif rounding_method == "probs":
            pooled_allele_counts = np.random.binomial(n=haploid_counts[i],
                                                      p=pop_pooled_allele_freqs,
                                                      size=num_of_polymorphisms)
        else:
            print("Error: Invalid rounding method.")
            sys.exit()

    # Construct numpy array containing SFS from pooled_allele_counts_by_population
    incremented_haploid_counts = [i + 1 for i in haploid_counts]
    sfs = np.zeros(incremented_haploid_counts)
    for i in range(num_of_polymorphisms):
        coords = []
        for pop_pooled_allele_counts in pooled_allele_counts_by_population:
            coords.append(pop_pooled_allele_counts[i])
        sfs[tuple(coords)] += 1

    return sfs

def main(vcf_file, popinfo_file, haploid_counts, poolseq_depth, rounding_method="counts", ploidy=2):
    sfs = get_pooled_folded_fs(vcf_file, popinfo_file, haploid_counts, poolseq_depth, rounding_method, ploidy):
    print(sfs)

# Set up working environment
working_dir = "/scratch/djb3ve/data/second_pipeline/"
vcf_file = working_dir + vcf_name
for subdir in ["vcfs/", "serialized_pooled_sfss/"]:
    if not os.path.exists(working_dir + subdir):
        os.makedirs(working_dir + subdir)

# Handle command line arguments
vcf_file = working_dir + "vcfs/" + sys.argv[1]
popinfo_name = working_dir + "popinfos/" + sys.argv[2]
haploid_counts = eval(sys.argv[3])
poolseq_depth = int(sys.argv[4])
rounding_method = sys.argv[5]
ploidy = int(sys.argv[6])

main()
