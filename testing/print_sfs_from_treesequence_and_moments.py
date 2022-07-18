import moments
import msprime
import tskit
import numpy as np

print("Done importing.")

mts = tskit.load("tps10.ts")
sfs_msprime = mts.allele_frequency_spectrum(sample_sets=[mts.samples()[:20], mts.samples()[20:]])

L = 1e6
for i in sfs_msprime:
    for j in i:
        print(int(j * L), end="\t")
    print()
print()

print(np.sum(sfs_msprime * L))

vcf_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/vcfs/two_pop_split_demography_n10_seed1.vcf"
popinfo_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/popinfos/two_pop_split_demography_n10_popinfo.txt"
data_dict = moments.Misc.make_data_dict_vcf(vcf_file, popinfo_file)
sfs_moments = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"], projections=[20, 20], polarized=False)

for i in sfs_moments:
    for j in i:
        print(j, end="\t")
    print()

print(np.sum(sfs_moments))
