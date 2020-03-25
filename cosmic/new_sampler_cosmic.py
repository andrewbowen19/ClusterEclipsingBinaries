"""Updating our cosmic_sampler function to include hs period cutoff, want to use only initial binaries"""

from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.sample.sampler import independent
from cosmic.sample.sampler import multidim
from cosmic.evolve import Evolve
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, \
'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, \
'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 3, 'ceflag': 0, 'eddfac': 1.0, 'merger': 0, 'ifflag': 0, \
'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'ppsn': 1,\
 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,\
  'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 0, \
  'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsnp' : 2.5, 'ecsn_mlow' : 1.6, 'aic' : 1, 'sigmadiv' :-20.0}


def cluster_bin_sampler(age, Nbin, Z, sigma, random_seed):
	"""Creates and evolves a set of binaries with given 
	age (to evolve to), number of binaries, metallicity, and velocity dispersion. 
	Will later loop through globular and open clusters and apply this for loop"""
	n_grid = Nbin

	# Initial (input) binares -- using sampler method from cosmic #1234 - random seed
	InitialBinaries, sampled_mass, n_sampled = InitialBinaryTable.sampler('multidim',\
	 [11], [11], random_seed, 1, 'delta_burst', age, Z, Nbin)#, porb_hi = [3])

	# Inclination and omega values
	inc = np.arccos(2.*np.random.uniform(0,1,Nbin) - 1.)
	omega = np.random.uniform(0,2*np.pi,Nbin)
	OMEGA = np.random.uniform(0,2*np.pi,Nbin)
	print(InitialBinaries)
# Making Input variables global (for plotting later)
	global p_i, m1_i, m2_i, ecc_i, tphysf_i, Z_i
	# Input binary params (for plotting later)
	p_i = InitialBinaries['porb']
	m1_i = InitialBinaries['mass1_binary']
	m2_i = InitialBinaries['mass2_binary']
	ecc_i = InitialBinaries['ecc']
	tphysf_i = InitialBinaries['tphysf']
	Z_i = InitialBinaries['metallicity']
	
	# Evolving input binaries
	bpp, bcm, initC  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

	# Making returned variables global
	global p_f, m1_f, m2_f, ecc_f, tphysf_f, r1_f, r2_f, lum1_f, lum2_f, Teff1, Teff2
	# Evolved Binary Params (made global for plotting later, can )
	p_f = bcm['porb']
	m1_f = bcm['mass_1']
	m2_f = bcm['mass_2']
	ecc_f = bcm['ecc']
	tphysf_f = bcm['tphys']
	r1_f = bcm['rad_1']
	r2_f = bcm['rad_2']
	lum1_f = bcm['lumin_1']
	lum2_f = bcm['lumin_2']
	Teff1 = bcm['teff_1']#Effective temperature - just in case
	Teff2 = bcm['teff_2']



	return (InitialBinaries, bcm)

bin_it = cluster_bin_sampler(13000, 1000, 0.02, 1, np.random.randint(1,100))
print(bin_it)













