import dadi

# Setup
vcf = "/scratch/djb3ve/data/test_poolseq.vcf"
popinfo = "/scratch/djb3ve/data/test_poolseq_popinfo.txt"

# Load VCF and popinfo file into "data_dict" in order to make SFS "fs"
data_dict = dadi.Misc.make_data_dict_vcf(vcf, popinfo)
fs = dadi.Spectrum.from_data_dict(data_dict,
                                  pop_ids=["pop1", "pop2", "pop3", "pop4"],
                                  projections=[3 for i in range(4)],
                                  polarized=False)

print(fs)
