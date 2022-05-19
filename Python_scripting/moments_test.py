import moments
print(-1)

input_vcf = "/scratch/djb3ve/data/small_vcf.vcf"
input_popinfo = "/scratch/djb3ve/data/popinfo.txt"

print(0)

data_dict = moments.Misc.make_data_dict_vcf(input_vcf, input_popinfo)

print(1)

fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"],
                                     projections=[5], polarized=True)

print (2)
