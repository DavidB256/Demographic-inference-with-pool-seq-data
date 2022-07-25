import msprime
import os
import yaml

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

if __name__ == "__main__":
    # Make necessary subdirectories "vcfs" and "popinfos" of "data_dir" if they do not already exist
    for subdirectory in ["vcfs/", "popinfos/"]:
        if not os.path.exists(yd["pipeline_params"]["data_dir"] + subdirectory):
            os.makedirs(yd["pipeline_params"]["data_dir"] + subdirectory)

    # Iterate through sample sizes, generating "num_of_replicates" many VCF files per sample size
    for sample_size in yd["pipeline_params"]["sample_sizes"]:
        print("Starting iteration with sample size %d." % sample_size)

        samples = {"pop0": sample_size, "pop1": sample_size}
        dem = msprime.Demography()
        dem.add_population(name="ancestral", initial_size=yd["dem_params"]["ancestral_population_size"])
        dem.add_population(name="pop0", initial_size=yd["dem_params"]["split_population_size"])
        dem.add_population(name="pop1", initial_size=yd["dem_params"]["split_population_size"])
        dem.add_population_split(time=yd["dem_params"]["t_split"], derived=["pop0", "pop1"], ancestral="ancestral")
        dem.set_symmetric_migration_rate(["pop0", "pop1"], yd["dem_params"]["mig_rate"])

        # Write popinfo file, which is used by moments and the pool-seq noising algorithm
        popinfo_file = f"{yd['pipeline_params']['data_dir']}popinfos/n{sample_size}_popinfo.txt"
        with open(popinfo_file, "w") as f:
            for i in range(sample_size * 2):
                f.write("tsk_%d\tpop%d\n" % (i, 0 if i < sample_size else 1))

        for rep in range(yd["pipeline_params"]["num_of_demography_sim_replicates"]):
            seed = rep + 1
            ts = msprime.sim_ancestry(samples=samples, demography=dem,
                                      sequence_length=yd["dem_params"]["seq_len"],
                                      recombination_rate=yd["dem_params"]["recombination_rate"],
                                      ploidy=yd["dem_params"]["ploidy"],
                                      random_seed=seed)
            mts = msprime.sim_mutations(ts, rate=yd["dem_params"]["mutation_rate"], random_seed=seed)

            # Serialize TreeSequence "mts"
            with open(f"{yd['pipeline_params']['data_dir']}tss/n{sample_size}_mseed{seed}.ts", "w") as f:
                mts.dump(f)
            # Export VCF created from TreeSequence "mts"
            with open(f"{yd['pipeline_params']['data_dir']}vcfs/n{sample_size}_mseed{seed}.vcf", "w") as f:
                mts.write_vcf(f)