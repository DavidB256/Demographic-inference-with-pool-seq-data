import msprime
import sys
import os.path
import os

class Demography_plus:
    def __init__(self, dem, dem_name, samples,
                 mu=2.9e-6, rho=1.25e-7, seq_len=1e5, pop_size=100,
                 ploidy=2, mig_rate=1e-2, t_split=10):
        # Set instance variables equal to constructor arguments
        self.dem = dem
        self.dem_name = dem_name
        self.samples = samples
        self.mu = mu
        self.rho = rho
        self.seq_len = seq_len
        self.pop_size = pop_size
        self.ploidy = ploidy
        self.mig_rate = mig_rate
        self.t_split = t_split

    # Returns tree sequence object from ancestry and mutation simulations performed
    # with demography and parameters specified by instance variables
    def get_ts_with_muts(self, random_seed=1):
        ts = msprime.sim_ancestry(samples=self.samples, demography=self.dem,
                                  sequence_length=self.seq_len,
                                  recombination_rate=self.rho,
                                  ploidy=self.ploidy, random_seed=random_seed)
        ts = msprime.sim_mutations(ts, rate=self.mu, random_seed=random_seed)
        return ts

    # Writes popinfo file into "output_dir" from "samples" instance variable
    def write_popinfo(self, output_dir):
        with open(output_dir + self.dem_name + "_popinfo.txt", "w") as f:
            sample_counter = 0
            # Samples is a dictionary whose keys are population names and whose values
            # are organism sample sizes
            for sample in self.samples:
                for i in range(self.samples[sample]):
                    f.write("tsk_%d\t%s\n" % (sample_counter, sample))
                    sample_counter += 1

    # Writes "iterations" many VCF files created from tree sequences created by
    # "get_ts_with_muts" with different random seeds. Writes the popinfo file if it
    # does not yet exist.
    def write_vcf_and_popinfo(self, output_dir, iterations=1):
        for i in range(iterations):
            with open(f"{output_dir}{self.dem_name}_seed{i}.vcf", "w") as f:
                self.get_ts_with_muts(random_seed=i).write_vcf(f)

        if not os.path.exists(output_dir + self.dem_name + "_popinfo.txt"):
            self.write_popinfo(output_dir)

    # Add lines to the pipeline instructions file to be read into poolseq_sfs_script.R
    # by its wrapper, one line per level of pool-seq depth
    def append_pipeline_instruction(self, instructions_output, poolseq_depths):
        # Format popinfo and haploid_counts as R vectors to be evaluated in the R script
        popinfo_r_vector = "c("
        haploid_counts_r_vector = "c("
        for sample in self.samples:
            popinfo_r_vector += str(self.samples[samples]) + ", "
            haploid_counts_r_vector += str(self.samples[samples] * self.ploidy) + ", "
        popinfo_r_vector = popinfo_r_vector[:-1] + ")"
        haploid_counts_r_vector = haploid_counts_r_vector[:-1] + ")"

        with open(instructions_output, "w+") as f:
            for poolseq_depth in poolseq_depths:
                f.write(f"{output_dir}{self.dem_name}_seed{i}.vcf" + "\t" + \
                        "\t" + popinfo_r_vector + "\t" + haploid_counts_r_vector + \
                        "\t" + str(poolseq_depth) + "\t" + instructions_output)

output_dir = "/scratch/djb3ve/data/first_models/"
instructions_output = output_dir + "pipeline_instructions.txt"
sample_sizes = [10, 50, 100, 200, 500]
poolseq_depths = [5, 10, 40, 70, 100]
iterations = 10

# Remove pre-existing instructions file and write header
with open(instructions_output, "w") as f:
    f.write("#VCF_file\tpopinfo\thaploid_counts\tpoolseq_depth\toutput_file")

for sample_size in sample_sizes:
    # control
    control_demography = msprime.Demography()
    control_demography.add_population(name="pop0", initial_size=100)
    dem_plus = Demography_plus(control_demography, "control_demography_n{sample_size}",
                               {"pop0": sample_size})
    dem_plus.write_vcf_and_popinfo(output_dir, iterations)
    dem_plus.append_pipeline_instruction(instructions_output, poolseq_depths)


    # two_pop_split
    two_pop_split_demography = msprime.Demography()
    two_pop_split_demography.add_population(name="pop0", initial_size=100)
    two_pop_split_demography.add_population(name="pop1", initial_size=100)
    two_pop_split_demography.add_population(name="ancestral", initial_size=200)
    two_pop_split_demography.add_population_split(time=t_split, derived=["pop0", "pop1"], ancestral="ancestral")
    two_pop_split_demography.set_symmetric_migration_rate(["pop0", "pop1"], 1e-2)
    dem_plus = Demography_plus(two_pop_split_demography, "two_pop_split_demography_n{sample_size}",
                               {"pop0": sample_size, "pop1": sample_size})
    dem_plus.write_vcf_and_popinfo(output_dir, iterations)
    dem_plus.append_pipeline_instruction(instructions_output, poolseq_depths)
