import moments

data_dict = moments.Misc.make_data_dict_vcf(vcf_filename, popinfo_filename)
fs = moments.Spectrum.from_data_dict(data_dict)
