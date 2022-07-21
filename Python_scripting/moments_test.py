import moments

# Setup
input_vcf = "/scratch/djb3ve/data/3island_small_model.vcf"
input_popinfo = "/scratch/djb3ve/data/popinfo_3island_small_model.txt"

# Load VCF and popinfo file into "data_dict" in order to make SFS "fs"
data_dict = moments.Misc.make_data_dict_vcf(input_vcf, input_popinfo)
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1", "pop2"],
                                     projections=[4, 4, 4], polarized=False)

for i in data_dict:
    print(i)
    print(data_dict[i])

print(fs)
