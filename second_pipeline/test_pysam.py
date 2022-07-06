import pandas as pd

vcf_df = pd.read_table("/scratch/djb3ve/data/second_pipeline/vcfs/two_pop_split_demography_n10_seed1.vcf", header=5)
popinfo_df = pd.read_table("/scratch/djb3ve/data/second_pipeline/popinfos/two_pop_split_demography_n10_popinfo.txt", header=None)

print(vcf_df[["tsk_0", "tsk_1"]])
print()

for i in popinfo_df[1].unique():
    x = popinfo_df[popinfo_df[1] == i][0]
    print(vcf_df[x])
    print()
