import msprime

rr_exps = [-6, -7, -8, -9, -10]
mu_exps = [-6, -5, -4, -3, -2]

for i, rr_exp in enumerate(rr_exps):
    for j, mu_exp in enumerate(mu_exps):
        ancestry = msprime.sim_ancestry(50, sequence_length=1_000_000, random_seed=1, recombination_rate=10 ** rr_exp)
        mut_ancestry = msprime.sim_mutations(ancestry, rate=10 ** mu_exp, random_seed=1)
        
        with open("%drr_%dmu_out.vcf" % (-rr_exp, -mu_exp), "w") as f:
            mut_ancestry.write_vcf(f)
        print("%d\t%d" % (i, j))
