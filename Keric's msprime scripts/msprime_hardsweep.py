import numpy as np
import msprime
import IPython
from IPython.display import SVG
import os
import math
import sys

#setting current working directory
os.chdir("/Users/kericlamb/Desktop/Documents/Research/Drosophila_phylogeography/msprime_neutral/")

#write vcf with given pre-determined scenario

#set up variables
popA_size = 15_000 #population size individuals are drawn from, held constant
sample_num=100 #individuals sampled
mu=2.8e-9 #mutation rate
rr=1.25e-7 #recombination rate
L=1e6 #length of genome
pl=2 #ploidy
sel=0.25 #selection coefficient of allele
dt=1e-6 #dt is the small increment of time for stepping through the sweep phase of the model. a good rule of thumb is for this to be approximately 1/40N or smaller.

#write demography
demography = msprime.Demography()
demography.add_population(name="A", initial_size=popA_size)

#set up sweep model
sweep_model = msprime.SweepGenicSelection(
    position=L / 2,  # middle of chrom
    start_frequency=1.0 / (2 * popA_size),
    end_frequency=1.0 - (1.0 / (2 * popA_size)),
    s=sel,
    dt=dt,
)

#choose sample size
sample = [msprime.SampleSet(sample_num, population= "A", ploidy=pl)] #randomly samples from population "A"

#simulate demographic history
ts = msprime.sim_ancestry(samples=sample, 
                          demography=demography, 
                          sequence_length=L, 
                          recombination_rate=rr, 
                          model=[sweep_model, "hudson"], 
                          ploidy=pl)
#simulate random mutations during history
mts = msprime.sim_mutations(ts, rate=mu)

#writes vcf
with open("/Users/kericlamb/Desktop/hardsweep.vcf", "w") as vcf_file:
                      mts.write_vcf(vcf_file)
