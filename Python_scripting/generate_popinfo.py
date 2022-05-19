import sys

# Handle command line arguments
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()
sample_size = int(sys.argv[1])

output_file = "/scratch/djb3ve/data/popinfo.txt"

with open(output_file, "w") as out:
    for i in range(sample_size):
        sample_name = "tsk_" + str(i)
        out.write("%s\tpop0" % (sample_name))
