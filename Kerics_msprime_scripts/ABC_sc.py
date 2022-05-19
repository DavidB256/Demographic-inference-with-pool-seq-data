import numpy as np
import msprime
import tskit
import os
import sys
import pandas as pd
import multiprocessing
import matplotlib.pyplot as plt
from IPython.display import SVG
import scipy
import seaborn as sns
from functools import reduce
import allel
import statistics

'''
required user input: 
number of simulations you want to run (sim_num) <- may expand in near future

'''

#defining constants
sim_num = int(sys.argv[0])
sim_range = range(sim_num) #number of simulations to run
k = 1000 #frequency of verbose checkups 

#setting up uniform distributions for pop size and time of split
#pop A
low_sizeA = 1e3
high_sizeA = 1e4
sample_A = 50
#pop B
low_sizeB = 1e3
high_sizeB = 1e4
sample_B = 50
#pop C
low_sizeC = 1e3
high_sizeC = 1e4
#Ts
low_Ts = 1e2
high_Ts = 1e4

#setting up exponential distribution center for migration
center_mig = 1e-4

#start secondary contact model section

#start summary statistic output file
PMmod=open('./ABC_model_sc.txt','w')
PMmod.write(
        str("sim_num")+'\t'+
        str("model")+'\t'+
        str("popA_size")+'\t'+
        str("popB_size")+'\t'+
        str("popC_size")+'\t'+
        str("T_split")+'\t'+
        str("T_mig")+'\t'+
        str("mig_rate")+'\t'+
        str("mean_fst")+'\t'+
        str("var_fst")+'\t'+
        str("mean_pairwise")+'\t'+
        str("var_pairwise")+'\t'+
        str("seq_div")+'\t'+
        str("taj_d1")+'\t'+
        str("taj_d2")+'\t'+
        str("diversity1")+'\t'+
        str("diversity2")+'\n')
PMmod.close()

for i in range(len(sim_range)):
    #model type
    model = str("split_sc")
    
    #defining run specific variables
    popA_size = np.random.uniform(low_sizeA, high_sizeA)
    popB_size = np.random.uniform(low_sizeB, high_sizeB)
    popC_size = np.random.uniform(low_sizeC, high_sizeC)
    T_split = np.random.uniform(low_Ts, high_Ts)
    T_mig = np.random.uniform(low_Ts, T_split) #setting up time of migration beginning for secondary contact
    mig_rate = np.random.exponential(center_mig)
    
    #create demography
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=popA_size) #generates descendant population A with size of 3k
    demography.add_population(name="B", initial_size=popB_size) #generates descendant population B with size of 2k
    demography.add_population(name="C", initial_size=popC_size) #generates ancestral population C with size of 5k
    demography.add_population_split(time=T_split, derived=["A", "B"], ancestral="C") #gives basic lineage orientation and time of split in generations
    demography.set_symmetric_migration_rate(["A", "B"], mig_rate)
    demography.add_symmetric_migration_rate_change(T_mig, ["A", "B"], 0)
    demography.sort_events()
    
    #set up sample draws from pop A and B
    sample_set_1 = msprime.SampleSet(sample_A, population= "A", ploidy=2) #sets ploidy and sample size
    sample_set_2 = msprime.SampleSet(sample_B, population= "B", ploidy=2) #sets ploidy and sample size 
    samples=[sample_set_1,sample_set_2] #concatenates

    #generates tree sequences
    ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=1e6, recombination_rate=2e-8, model="hudson", ploidy=2) #generates basic tree sequence with known seq length, recomb rate, and demography
    mts = msprime.sim_mutations(ts, rate=2.8e-9) #incorporates random mutations

    #writes VCF
    with open("msprime_sim.vcf", "w") as vcf_file:
                          mts.write_vcf(vcf_file)

    #calculate summary statistics from vcf using scikit-allel
    
    #read in vcf generated above
    #read in meta dat (generalized to work for all sims-- only need to generate once)
    callset = allel.read_vcf('msprime_sim.vcf')
    meta = pd.read_csv("sim_vcf_meta.txt", sep='\t')
    
    #get positions and genotypes
    pos_all = allel.SortedIndex(callset['variants/POS']) #call positions
    genotype_all = allel.GenotypeChunkedArray(callset['calldata/GT']) #call genotypes and make GT array
    
    #set up populations
    pop1 = "pop1"
    pop2 = "pop2"
    n_samples_pop1 = np.count_nonzero(meta.population == pop1)
    n_samples_pop2 = np.count_nonzero(meta.population == pop2)
    
    subpops = {pop1: meta[meta.population == pop1].index, pop2: meta[meta.population == pop2].index, }
    
    #get allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    
    #sets up allele count arrays
    acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
    flt = acu.is_segregating() & (acu.max_allele() == 1) #grabs seg sites and gets rid of multivariant (>1 variant) sites for simplicity
    
    #sets up final positions (compressed for memory sake)
    pos = pos_all.compress(flt)
    #sets up individual allele count arrays for each population (necessary for calculating sum stats)
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
    #sets up final compressed genotype array
    genotype = genotype_all.compress(flt, axis=0)
    
    #now ready to calculate summary statistics
    
    #Fst first
    # sample indices for population 1
    pop1_idx = subpops[pop1]
    # sample indices for population 2
    pop2_idx = subpops[pop2]

    #estimates theta/Fst for each variant and each allele individually
    a, b, c = allel.weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], max_allele=1)
    snp_fst_wc = (a / (a + b + c))[:, 0]
    mean_fst = statistics.mean(snp_fst_wc)
    var_fst = statistics.mean(snp_fst_wc)
    
    #now mean pairwise differences between pops
    mean_pairwise_diff = allel.mean_pairwise_difference_between(ac1, ac2)
    mean_pairwise = statistics.mean(mean_pairwise_diff)
    var_pairwise = statistics.variance(mean_pairwise_diff)
    
    #sequence divergence
    seq_div = allel.sequence_divergence(pos, ac1, ac2)
    
    #Tajima's D
    taj_d1 = allel.tajima_d(ac1)
    taj_d2 = allel.tajima_d(ac2)
    
    #sequence diversity
    diversity1 = allel.sequence_diversity(pos, ac1)
    diversity2 = allel.sequence_diversity(pos, ac2)

    #prints update on how far along simulations are
    #prints if remainder equals 0, so pick a big number like 100 or 250 unless you want to drown in verbose updates
    if sim_range[i] % k == 0:
        print("finishing simulation number:", sim_range[i])
    
    PMmod=open('./ABC_model_sc.txt','a')
    PMmod.write(
            str(sim_range[i])+'\t'+
            str(model)+'\t'+
            str(popA_size)+'\t'+
            str(popB_size)+'\t'+
            str(popC_size)+'\t'+
            str(T_split)+'\t'+
            str(T_mig)+'\t'+
            str(mig_rate)+'\t'+
            str(mean_fst)+'\t'+
            str(var_fst)+'\t'+
            str(mean_pairwise)+'\t'+
            str(var_pairwise)+'\t'+
            str(seq_div)+'\t'+
            str(taj_d1)+'\t'+
            str(taj_d2)+'\t'+
            str(diversity1)+'\t'+
            str(diversity2)+'\n')
    PMmod.close()

print("done simulating split with secondary contact scenario")