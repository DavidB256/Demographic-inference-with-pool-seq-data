import sys

# Handle command line arguments
if len(sys.argv) < 2:
    print("Error: One command line argument is required.")
    sys.exit()
input_vcf = sys.argv[1]

out_dir = "/scratch/djb3ve/data/"

with open(input_vcf, "r") as vcf:
    vcf_file_name_wo_ext = input_vcf.split(".")[-2].split("/")[-1]
    with open(out_dir + vcf_file_name_wo_ext + "_popinfo.txt", "w") as out:
        out.write("SAMPLE_NAME\tPOP_NAME\n")
        line = vcf.readline()
        while line:
            if line[0] != "#":
                pos = line.split()[1]
                out.write("%s\t0\n" % (pos))
            line = vcf.readline()
