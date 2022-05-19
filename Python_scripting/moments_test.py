import moments
print(-1)
input_vcf = "/scratch/djb3ve/data/3islands_4mig.vcf"

print(0)
data_dict = moments.Misc.make_data_dict_vcf(input_vcf)

print(1)

fs = moments.Spectrum.from_data_dict(data_dict)

print (2)
