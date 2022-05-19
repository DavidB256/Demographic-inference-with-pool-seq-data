import dadi

input_vcf = "/scratch/djb3ve/data/small_vcf.vcf"
input_popinfo = "/scratch/djb3ve/data/popinfo.txt"
data_dict = dadi.Misc.make_data_dict_vcf(input_vcf, input_popinfo)
fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"],
                                     projections=[5], polarized=False)

pi = fs.pi()

print(pi)
