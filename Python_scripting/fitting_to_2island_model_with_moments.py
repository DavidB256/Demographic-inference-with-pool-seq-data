import moments
import dadi

# Import VCF file from prior msprime simulation and popinfo file
vcf = "/scratch/djb3ve/data/2island_1mig_model.vcf"
popinfo = "/scratch/djb3ve/data/popinfo_file_for_2island_model_10n.txt"
output = "/scratch/djb3ve/yuh.txt"
pts = 100
ns = [10, 10]

print("Setup complete.")

data_dict = dadi.Misc.make_data_dict_vcf(vcf, popinfo)

print("VCF imported.")

fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids=["pop0", "pop1"],
                                     projections=[5], polarized=False)

print("VCF converted to SFS.")

def two_island_admixture(params, ns):
    nu1, nu2, T, m = params
    xx = moments.default_grid(pts)

    phi = moments.PhiManip.phi_1D(xx)
    phi = moments.PhiManip.phi_1D_to_2D(xx, phi)
    phi = moments.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)
    fs = moments.Spectrum.from_phi(phi, ns, (xx, xx))

    return fs

print("Model function defined.")

mig_rates = [i * 0.01 for i in range(25)]
out_f = open(output, "w")

for mig_rate in mig_rates:
    params = [100, 100, 0, mig_rate]
    popt = moments.Inference.optimate_log(params, fs, two_island_admixture, pts)
    model = two_island_admixture(popt, ns)

    ll_model = moments.Inference.ll_multinom(model, fs)
    print("%f\t%f" % (mig_rate, ll_model))
    out_f.write(print("%f\t%f" % (mig_rate, ll_model)))

out_f.close()
