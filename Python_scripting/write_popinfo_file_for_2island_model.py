import sys

out_dir = "/scratch/djb3ve/data/"

# Handle command line arguments
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()
sample_size = int(sys.argv[1])

with open(out_dir + "popinfo_file_for_2island_model_%dn.txt"
          % (sample_size), "w") as out:
    for i in range(sample_size):
        out.write("tsk_%d\tpop0\n" % (i))
    for i in range(sample_size):
        out.write("tsk_%d\tpop1\n" % (i + sample_size))
