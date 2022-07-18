import numpy as np
import abcpy
from abcpy.continuousmodels import Normal as Gaussian, Uniform
from abcpy.statistics import Identity
from abcpy.distances import LogReg
from abcpy.perturbationkernel import DefaultKernel
from abcpy.backends import BackendDummy as Backend
from abcpy.inferences import PMCABC

height_obs = [160.82499176, 167.24266737, 185.71695756, 153.7045709, 163.40568812, 140.70658699, 169.59102084,
              172.81041696, 187.38782738, 179.66358934, 176.63417241, 189.16082803, 181.98288443, 170.18565017,
              183.78493886, 166.58387299, 161.9521899, 155.69213073, 156.17867343, 144.51580379, 170.29847515,
              197.96767899, 153.36646527, 162.22710198, 158.70012047, 178.53470703, 170.77697743, 164.31392633,
              165.88595994, 177.38083686, 146.67058471763457, 179.41946565658628, 238.02751620619537,
              206.22458790620766, 220.89530574344568, 221.04082532837026, 142.25301427453394, 261.37656571434275,
              171.63761180867033, 210.28121820385866, 237.29130237612236, 175.75558340169619, 224.54340549862235,
              197.42448680731226, 165.88273684581381, 166.55094082844519, 229.54308602661584, 222.99844054358519,
              185.30223966014586, 152.69149367593846, 206.94372818527413, 256.35498655339154, 165.43140916577741,
              250.19273595481803, 148.87781549665536, 223.05547559193792, 230.03418198709608, 146.13611923127021,
              138.24716809523139, 179.26755740864527, 141.21704876815426, 170.89587081800852, 222.96391329259626,
              188.27229523693822, 202.67075179617672, 211.75963110985992, 217.45423324370509]
mu = Uniform([[150], [200]], name="mu")
sigma = Uniform([[5], [25]], name="sigma")
height = Gaussian([mu, sigma], name='height')

statistics_calculator = Identity(degree=2, cross=False)
distance_calculator = LogReg(statistics_calculator, seed=42)

kernel = DefaultKernel([mu, sigma])
backend = Backend()
sampler = PMCABC([height], [distance_calculator], backend, kernel, seed=1)

steps = 10
eps_arr = np.array([.75])
n_sample = 10
n_samples_per_param = 10
epsilon_percentile = 10

journal = sampler.sample([height_obs], steps, eps_arr, n_sample, n_samples_per_param, epsilon_percentile)
print(journal.get_parameters())
print(journal.get_weights())
