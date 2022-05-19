import dadi

print(-1)

input_vcf = "/scratch/djb3ve/data/small_vcf.vcf"
input_popinfo = "/scratch/djb3ve/data/small_vcf_popinfo.txt"

print(0)

data_dict = dadi.Misc.make_data_dict_vcf(input_vcf, input_popinfo)

print(1)

fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["0"],
                                     projections=[5], polarized=False)

print (2)
