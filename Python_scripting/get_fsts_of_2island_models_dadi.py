import dadi

vcfs = ["/scratch/djb3ve/data/2island_%dmig_model.vcf" % (i + 1)
        for i in range(4)]

# need to generate multiple popinfo files, but they can all be identical????

input_popinfo = "/scratch/djb3ve/data/popinfo.txt"
data_dict = dadi.Misc.make_data_dict_vcf(input_vcf, input_popinfo)
fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0"],
                                     projections=[5], polarized=False)

pi = fs.pi()

print(pi)
