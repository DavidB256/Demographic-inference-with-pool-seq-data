import sys
import os
import numpy as np
import yaml
import re
import msprime
import tskit

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)


def add_noise_to_sfs(sfs, poolseq_depth, method):
    haploid_counts = [i - 1 for i in sfs.shape]
    noised_sfs = np.zeros(sfs.shape)
    it  = np.nditer(sfs, flags=['multi_index'])

    for i in it:
        freqs = [index / haploid_counts[j] for j, index in enumerate(it.multi_index)]
        for allele in range(i):
            # Generate an allele's frequencies in each subpopulation.
            # Uses same method as "sample.alleles(mode="coverage")" in Thomas Taus' poolSeq R pakage.
            # (https://github.com/ThomasTaus/poolSeq/blob/master/R/simaf.R).
            # Lacks support for specifying pool-seq coverage at each site with vector argument for "size" parameter.
            # Lacks support for mode="individuals".
            coverages = np.random.poisson(lam=poolseq_depth, size=sfs.ndim)
            coverages = np.where(coverages == 0, 1, coverages)
            pooled_freqs = np.random.binomial(n=coverages, p=freqs, size=sfs.ndim) / coverages

            # Convert frequencies to counts via a rounding method.
            # The rounding methods are taken from genomalicious' "dadi_inputs_pools" function
            # (https://rdrr.io/github/j-a-thia/genomalicious/man/dadi_inputs_pools.html).
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

if __name__ == "__main__":
    # Handle command line arguments
    ts_name = sys.argv[1]
    ts_file = yd["pipeline_params"]["data_dir"] + "tss/" + sys.argv[1]

    sample_size = int(re.findall(r'\d+', ts_name)[0])

    # Use msprime to construct SFS from VCF.
    ts = tskit.load(ts_file)
    sfs = ts.allele_frequency_spectrum(sample_sets=[ts.samples()[:sample_size * yd["dem_params"]["ploidy"]],
                                                    ts.samples()[sample_size * yd["dem_params"]["ploidy"]:]])
    sfs *= yd["dem_params"]["seq_len"]
    sfs = sfs.astype(int)

    # Serialize non-noised SFS
    output_file = yd["pipeline_params"]["output_dir"] + ts_name[:-3] + "_depth0_pseed1"
    np.save(output_file, sfs)

    for poolseq_depth in yd["pipeline_params"]["poolseq_depths"]:
        # Iterate once per replicate, incrementing the seed between replicates
        for seed in range(1, yd["pipeline_params"]["num_of_demography_sim_replicates"]+1):
            np.random.seed(seed)

            # Determine SFS
            noised_sfs = add_noise_to_sfs(sfs, poolseq_depth, yd["pipeline_params"]["rounding_method"])
            # The file extension is not included in "output_file" because "np.save" automatically adds the ".npy" extension.
            output_file = yd["pipeline_params"]["output_dir"] + ts_name[:-3] + "_depth" + str(poolseq_depth) + "_pseed" + str(seed)
            # Serialize SFS
            np.save(output_file, noised_sfs)





