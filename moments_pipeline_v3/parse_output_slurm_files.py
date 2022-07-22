import os
import yaml

# Collects moments output data from Slurm output files into "output_file".

# Import config YAML file into global dictionary
with open("config.yaml", "r") as f:
    yd = yaml.safe_load(f)

data_dir = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/moments_pipeline_v3/output_slurm/"
output_file = yd["pipeline_params"]["output_file"]

with open(output_file, "a+") as g:
    for file in os.listdir(data_dir):
        with open(file, "r") as f:
            line = f.readline()
        g.write(line)
