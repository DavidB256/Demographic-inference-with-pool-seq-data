import msprime
import sys
import os

# Given the parameters at the bottom of this script and the demography model described
# in main(), this script writes VCF files replicates of VCF files at different sample sample sizes
# and with different random seeds with msprime and writes corresponding popinfo files
# into the "second_pipeline" directory.

# Objects of this class have demography objects and metadata for simulation and
# the pipeline as instance variables
class Demography_plus:
    # Values for mu and rho are those of D. melanogaster
    def __init__(self, dem, dem_name, samples,
                 mu, rho, seq_len, ploidy=2):
        # Set instance variables equal to constructor arguments
        self.dem = dem
        self.dem_name = dem_name
        self.samples = samples
        self.mu = mu
        self.rho = rho
        self.seq_len = seq_len
        self.ploidy = ploidy

    # Returns tree sequence object from ancestry and mutation simulations performed
    # with demography and parameters specified by instance variables
    def get_ts_with_muts(self, random_seed=1):
        ts = msprime.sim_ancestry(samples=self.samples, demography=self.dem,
                                  sequence_length=self.seq_len,
                                  recombination_rate=self.rho,
                                  ploidy=self.ploidy, random_seed=random_seed)
        return msprime.sim_mutations(ts, rate=self.mu, random_seed=random_seed)

    # Writes popinfo file into "output_dir" from "samples" instance variable
    def write_popinfo(self, output_file):
        with open(output_file, "w") as f:
            sample_counter = 0
            # Samples is a dictionary whose keys are population names and whose values
            # are organism sample sizes
            for sample in self.samples:
                for i in range(self.samples[sample]):
                    f.write("tsk_%d\t%s\n" % (sample_counter, sample))
                    sample_counter += 1

    # Writes "num_of_replicates" many VCF files created from tree sequences created by
    # "get_ts_with_muts" with different random seeds. Writes the popinfo file if it
    # does not yet exist.
    def write_vcf_and_popinfo(self, num_of_replicates, output_dir):
        # Need to iterate with 1-based counting for seeds because a seed of 0 is invalid
        for seed in range(1, num_of_replicates+1):
            with open(f"{output_dir}vcfs/{self.dem_name}_seed{seed}.vcf", "w") as f:
                self.get_ts_with_muts(random_seed=seed).write_vcf(f)

        if not os.path.exists(f"{output_dir}popinfos/{self.dem_name}_popinfo.txt"):
            self.write_popinfo(f"{output_dir}popinfos/{self.dem_name}_popinfo.txt")

def main(dem_params, output_dir, sample_sizes, num_of_replicates):
    # Make necessary subdirectories "vcfs" and "popinfos" of "output_dir" if they do not already exists
    for subdirectory in ["vcfs", "popinfos"]:
        if not os.path.exists(output_dir + subdirectory):
            os.makedirs(output_dir + subdirectory)

    # Iterate through sample sizes, generating "num_of_replicates" many VCF files per sample size
    for sample_size in sample_sizes:
        print("Starting iteration with sample size %d." % sample_size)
        # Create Demography_plus object for two_pop_split model
        two_pop_split_demography = msprime.Demography()
        two_pop_split_demography.add_population(name="ancestral",
                                                initial_size=dem_params["ancestral_population_size"])
        two_pop_split_demography.add_population(name="pop0",
                                                initial_size=dem_params["split_population_size"])
        two_pop_split_demography.add_population(name="pop1",
                                                initial_size=dem_params["split_population_size"])
        two_pop_split_demography.add_population_split(time=dem_params["t_split"],
                                                      derived=["pop0", "pop1"], ancestral="ancestral")
        two_pop_split_demography.set_symmetric_migration_rate(["pop0", "pop1"], dem_params["mig_rate"])
        dem_plus = Demography_plus(two_pop_split_demography, f"two_pop_split_demography_n{sample_size}",
                                   {"pop0": sample_size, "pop1": sample_size},
                                   dem_params["mu"], dem_params["rho"], dem_params["seq_len"], dem_params["ploidy"])
        dem_plus.write_vcf_and_popinfo(num_of_replicates, output_dir)

# Demography parameters
dem_params = {"ancestral_population_size": 2e5,
              "split_population_size": 1e5,
              "t_split": 10,
              "mig_rate": 1e-2,
              "mu": 2.9e-9,
              "rho": 1.25e-7,
              "seq_len": 1e6,
              "ploidy": 2}

# Procedural variables
output_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/"
sample_sizes = [10, 50, 100, 200, 400]
num_of_replicates = 3

main(dem_params, output_dir, sample_sizes, num_of_replicates)
