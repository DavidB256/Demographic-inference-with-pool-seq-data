import sys
import os
import numpy as np
import pandas as pd

working_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"
vcf_file = working_dir + "vcfs/two_pop_split_demography_n10_seed1.vcf"
popinfo_file = working_dir + "popinfos/two_pop_split_demography_n10_popinfo.txt"

haploid_counts = [20, 20]

def get_allele_count_from_vcf_row(row):
    cumulative_count = 0
    for polymorphism in row:
        allele_list = polymorphism.split("|")
        alt_allele_count = len(allele_list) - allele_list.count("0")
        cumulative_count += alt_allele_count
    return cumulative_count

vcf_df = pd.read_table(vcf_file, header=5)
# Remove non-biallelic polymorphisms
vcf_df = vcf_df[vcf_df.ALT.apply(lambda x: len(x) == 1)]
vcf_df.reset_index(drop=True, inplace=True)
popinfo_df = pd.read_table(popinfo_file, header=None)

# Number of rows of "vcf_df" is the number of polymorphic sites recorded in the VCF
num_of_polymorphisms = vcf_df.shape[0]

# Subdivide columns of VCF file by population, with one population per element of
# "vcf_dfs_by_population"
vcf_dfs_by_population = []
for pop_id in popinfo_df[1].unique():
    tsks = popinfo_df[popinfo_df[1] == pop_id][0]
    vcf_dfs_by_population.append(vcf_df[tsks])

# Obtain allele counts from subsets of rows of VCF file corresponding to each population
allele_counts_by_population = pd.DataFrame()
for i, pop_vcf_df in enumerate(vcf_dfs_by_population):
    allele_counts = pop_vcf_df.apply(get_allele_count_from_vcf_row, axis=1)
    allele_counts_by_population[i] = allele_counts

total_haploids_sampled = np.sum(haploid_counts)
def fold_allele_counts(row):
    if np.sum(row) > total_haploids_sampled / 2:
        for i, count in enumerate(row):
            row[i] = min(count, haploid_counts[i] - count)
    return row
folded_allele_counts_by_population = allele_counts_by_population.apply(fold_allele_counts, axis=1)

# Homebrew algorithm for converting allele counts into SFS, works for nD
incremented_haploid_counts = [i + 1 for i in haploid_counts]
sfs = np.zeros(incremented_haploid_counts)
for i, row in folded_allele_counts_by_population.iterrows():
    sfs[tuple(row)] += 1

# Mask corners, which represent alleles that are absent and ubiquitous, respectively
absent_allele_coords = [0 for i in haploid_counts]
sfs[tuple(absent_allele_coords)] = 0
sfs[tuple(haploid_counts)] = 0


# Replicate of scikit-allel's joint_sfs(), only works for 2D
"""
tmp = (allele_counts_by_population[0] * haploid_counts[1] + allele_counts_by_population[1])
sfs = np.bincount(tmp)
# Mask corners
sfs[0] = 0
sfs[-1] = 0
# Reshape from 1D to nD
sfs.resize(haploid_counts[0], haploid_counts[1])
"""

# Print SFS
for i in sfs:
    for j in i:
        print(j, end="\t")
    print()
print(np.sum(sfs))
