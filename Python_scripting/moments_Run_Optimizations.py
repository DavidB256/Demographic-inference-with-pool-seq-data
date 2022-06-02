'''
Usage: python moments_Run_Optimizations.py
This is meant to be a general use script to run moments to fit any model on an
afs/jsfs with one to five populations. The user will have to edit information about
their allele frequency spectrum and provide a custom model. The instructions are
annotated below, with a #************** marking sections that will have to be edited.
Several examples of how to use various arguments to control optimizations are shown.
This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary.
If you'd like to use this script for larger sets of models already available, please
look on the other repositories to see how to import models from external model scripts.
General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:

 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688
Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject.
 The chi-squared test statistic is calculated assuming the sfs is folded. If this is not
 true for your data set, this number will not be accurate. This could be edited in the
 'collect_results' function in the Optimize_Functions.py script for an unfolded spectrum.
-------------------------
Written for Python 2.7 and 3.7
Python modules required:
  -Numpy
  -Scipy
  -moments
-------------------------
Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated August 2019
'''

import sys
import os
import numpy
import moments
import pylab
from datetime import datetime
import Optimize_Functions

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************
snps = "/Users/portik/Dropbox/MOMENTS/moments_pipeline/Example_Data/moments_2pops_North_South_snps.txt"

#Create python dictionary from snps file
dd = moments.Misc.make_data_dict(snps)

#**************
#pop_ids is a list which should match the populations headers of your SNPs file columns
pop_ids=["North", "South"]

#**************
#projection sizes, in ALLELES not individuals
proj = [16, 32]

#Convert this dictionary into folded AFS object
#[polarized = False] creates folded spectrum object
fs = moments.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(proj))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#================================================================================
# Here is an example of using a custom model within this script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script:

 Optimize_Routine(fs, outfile, model_name, func, rounds, param_number, fs_folded=True,
                       reps=None, maxiters=None, folds=None, in_params=None,
                       in_upper=None, in_lower=None, param_labels=" "):
   Mandatory Arguments =
    fs:  spectrum object name
    outfile:  prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or
          from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).
   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order
'''


# Let's start by defining our model
def sym_mig(params, ns):
    """
    Split into two populations, with symmetric migration.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations)
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    sym_mig = numpy.array([[0,m], [m,0]])
    fs.integrate([nu1,nu2], T, m=sym_mig, dt_fac=0.01)
    return fs



'''
Example 1. Now let's use the function to run an optimization routine for our data and this model.
We need to specify the first six arguments in this function, but there are other options
we can also use if we wanted more control over the optimization scheme. We'll start with
the basic version here. The argument explanations are above. This would perform three
rounds of optimizations, using a default number of replicates for each round.
'''
#create a prefix to label the output files
prefix = "V1"
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True)




'''
Example 2. It is a good idea to include the labels of the parameters so they can get written to the
output file, otherwise you'll have to go back to the model each time you wanted to see their
order. The optional arguments require using the = sign to assign a variable or value to the argument.
'''
prefix = "V2"

plabels = "nu1, nu2, m, T"

Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True,
                                        param_labels = plabels)




'''
Example 3. Here is the same example but also including your own custom parameter bounds. Notice the
optional arguments can be placed in any order following the mandatory arguments.
'''
prefix = "V3"
plabels = "nu1, nu2, m, T"

upper = [20,20,10,15]
lower = [0.01,0.01,0.01,0.1]

Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True,
                                        param_labels = plabels, in_upper = upper, in_lower = lower)





'''
Example 4. You can also be very explicit about the optimization routine, controlling what happens
across each round. Let's keep the three rounds, but change the number of replicates,
the maxiter argument, and fold argument each time. We'll need to create a list of values
for each of these, that has three values within (to match three rounds).
'''
prefix = "V4"
plabels = "nu1, nu2, m, T"
upper = [20,20,10,15]
lower = [0.01,0.01,0.01,0.1]

reps = [10,20,50]
maxiters = [5,10,20]
folds = [3,2,1]
'''
Using these arguments will cause round one to have 10 replicates, use 3-fold perturbed
starting parameters, and a maxiter of 5 for the optimization algorithm steps. Round two
will have 20 replicates, use 2-fold perturbed starting parameters, and a maxiter of 10
for the optimization algorithm steps, and etc. for round three.
'''
Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True,
                                        param_labels = plabels, in_upper=upper, in_lower=lower,
                                        reps = reps, maxiters = maxiters, folds = folds)




'''
Example 5. It's also good run the optimization routine multiple times. Let's write a short
loop to do the above optimization routine five times. We will name the prefix based
on which point we are at, and include it within the loops. Note that when you use
the range argument in python it will go up to, but not include, the final number.
That's why I have written a range of 1-6 to perform this 5 times.
'''
plabels = "nu1, nu2, m, T"
upper = [20,20,10,15]
lower = [0.01,0.01,0.01,0.1]
reps = [10,20,50]
maxiters = [5,10,20]
folds = [3,2,1]

for i in range(1,6):
    prefix = "V5_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(fs, prefix, "sym_mig", sym_mig, 3, 4, fs_folded=True,
                                            param_labels = plabels, in_upper=upper, in_lower=lower,
                                            reps = reps, maxiters = maxiters, folds = folds)
