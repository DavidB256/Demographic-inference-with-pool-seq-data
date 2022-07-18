import numpy as np
npy_file = "/scratch/djb3ve/Demographic-inference-with-pool-seq-data/data/second_pipeline/sfss/two_pop_split_demography_n10_seed1_depth100_seed1.npy"

sfs = np.load(npy_file)

for i in sfs:
    for j in i:
        print(int(j), end="\t")
    print()
