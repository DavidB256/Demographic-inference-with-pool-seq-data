import moments

# Setup
input_vcf = "/scratch/djb3ve/data/small_vcf.vcf"
input_popinfo = "/scratch/djb3ve/data/popinfo_for_small_vcf.txt"

# Load VCF and popinfo file into "data_dict" in order to make SFS "fs"
data_dict = moments.Misc.make_data_dict_vcf(input_vcf, input_popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"],
                                     projections=[5], polarized=False)

# Calculate and print pi
pi = fs.pi()
print(fs)
