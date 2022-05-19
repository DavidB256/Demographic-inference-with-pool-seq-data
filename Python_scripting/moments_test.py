import moments

input_vcf = "/scratch/djb3ve/data/3islands_4mig.vcf"

data_dict = moments.Misc.make_data_dict_vcf(input_vcf)
fs = moments.Spectrum.from_data_dict(data_dict)
