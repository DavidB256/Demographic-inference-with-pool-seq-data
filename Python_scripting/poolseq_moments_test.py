import msprime
import numpy as np

input_fs = "/scratch/djb3ve/data/test_poolseq_sfs.txt"

with open(input_fs, "r") as f_fs:
    dims = []
    values = []
    line = f_fs.readline()
    while line[:-1].isnumeric():
        dims.append(int(line[:-1]))
        line = f_fs.readline()
    line = f_fs.readline()
    while line:
        values.append(int(line[:-1]))
        line = f_fs.readline()

arr = np.transpose(np.array(values).reshape(dims))
