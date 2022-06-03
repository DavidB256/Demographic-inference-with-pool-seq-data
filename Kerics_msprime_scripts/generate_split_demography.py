import numpy as np
import msprime
import os
import math
import os
import sys
import numpy as np

#write vcf with given pre-determined scenario

#set up variables
popA_size = 3_000
popB_size = 2_000
popC_size = 5_000
time_split = 1000
sample_pop1 = 50
sample_pop2 = 50

#write demography
demography = msprime.Demography()
demography.add_population(name="A", initial_size=popA_size)
demography.add_population(name="B", initial_size=popB_size)
demography.add_population(name="C", initial_size=popC_size)
demography.add_population_split(time=time_split, derived=["A", "B"], ancestral="C") #this makes A and B cleave off from C at a predetermined timepoint
demography.set_symmetric_migration_rate(["A", "B"], 3e-5)

#choose sample size
sample_set_1 = msprime.SampleSet(sample_pop1, population= "A", ploidy=2)
sample_set_2 = msprime.SampleSet(sample_pop2, population= "B", ploidy=2)
samples=[sample_set_1,sample_set_2]

#simulate demographic history
ts = msprime.sim_ancestry(samples=samples, demography=demography, sequence_length=1e6,
                          recombination_rate=2e-8, model="hudson", ploidy=2)

#simulate random mutations during history

mts = msprime.sim_mutations(ts, rate=1e-9)

#writes vcf
with open("/scratch/djb3ve/data/msprime_null_split.vcf", "w") as vcf_file:
                      mts.write_vcf(vcf_file)
