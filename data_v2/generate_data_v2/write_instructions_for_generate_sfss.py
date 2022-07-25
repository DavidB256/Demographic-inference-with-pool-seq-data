import yaml

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

# Write instructions file with columns:
# vcf_name, popinfo_name, poolseq_depth
with open("instructions_for_generate_sfss.txt", "w") as f:
    for sample_size in yd["pipeline_params"]["sample_sizes"]:
        popinfo_name = f"n{sample_size}_popinfo.txt"
        for seed in range(1, yd["pipeline_params"]["num_of_demography_sim_replicates"]+1):
            vcf_name = f"n{sample_size}_mseed{seed}.vcf"
            f.write("%s\t%s\n" % (vcf_name, popinfo_name))
