import msprime
import numpy as np
import sys

input_fs = "/scratch/djb3ve/data/" + sys.argv[1]

with open(input_fs, "r") as f_fs:
    dims = []
    line = f_fs.readline()
    while line[:-1].isnumeric():
        dims.append(int(line[:-1]))
        line = f_fs.readline()
    values = []
    line = f_fs.readline()
    while line:
        values.append(int(line[:-1]))
        line = f_fs.readline()

arr = np.transpose(np.array(values).reshape(dims))
print(arr)
