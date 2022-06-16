import moments
import numpy as np
import sys
import os
import re
import statistics

def two_pop_split(params, ns):
    return moments.Demographics2D.split_mig(params, ns)
