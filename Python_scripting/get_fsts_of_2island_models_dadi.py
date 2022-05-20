import dadi

data_dir = "/scratch/djb3ve/data/"

vcfs = ["%s2island_%dmig_model.vcf" % (data_dir, i + 1)
        for i in range(4)]
popinfo = data_dir + "popinfo_file_for_2island_model_10n.txt"

for vcf in vcfs:
    data_dict = dadi.Misc.make_data_dict_vcf(vcf, popinfo)
    fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                      projections=[10, 10], polarized=False)

    pi = fs.pi()
    print(pi)
