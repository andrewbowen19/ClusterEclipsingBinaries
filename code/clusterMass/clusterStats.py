'''
Testing script to produce stats on cluster params
Prints avg and summed values for each cluster type to screen
Also writes options to .txt files in output foler (located in current dir)
Could also be nice to have user input values desired
'''

import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from itertools import permutations

# Print values to screen (for OCs and GCs)
# Writes values to txt file as well
# Maybe include plots??? idk just messing around
# If some one passes 'OC' or 'GC' command line arg then maybe 

class clusterStats(object):
	'''
	Object that includes all the functionality in producing/printing cluster stats
	Generates output (.txt) files as well as plots
	clusterParams - list-like of cluster parameter keys desired 
	*Note that param keys are not the same as the column labels (see self.paramDict)
	clusterType: can be
	'''

	def __init__(self, clusterParams, clusterType):
		self.clusterParams = clusterParams
		self.clusterType = clusterType
		self.paramDict = {'dist':'dist[pc]', 'rhm':'rhm[pc]', 'mass':'mass[Msun]',
						  'age':'age[Myr]', 'z':'[Fe/H]', 'sigma':'sigma_v0_z[km/s]'}
		# Dictionary containing units for each param key
		self.unitsDict = {'dist':'pc', 'rhm':'pc', 'mass':'M_sol',
						  'age':'Myr', 'z':'Fe/H', 'sigma':'km/s'}
		# self.df = None


	def getClusterData(self):
		'''
		Reading in cluster data csv files
		need to be located in data sub-directory of working dir
		'''
		
		self.gc = pd.read_csv(os.path.join('data', 'GCdataForEBLSST.csv'), header=0)
		self.oc = pd.read_csv(os.path.join('data', 'OCdataForEBLSST.csv'), header=0)

		# Printing data to screen - uncomment if desired
		# print('Globular Cluster Data: \n', self.gc)
		# print('Open Cluster Data: \n', self.oc)

		# Setting up 
		if 'gc' in self.clusterType.lower() or 'globular' in self.clusterType.lower():
			self.df = self.gc
			self.cluster_dfs = [self.gc]

		if 'oc' in self.clusterType.lower() or 'open' in self.clusterType.lower():
			self.df = self.oc
			self.cluster_dfs = [self.oc]

		# If clusterType is passed as 'All' or 'Both' or combined
		if 'all' in self.clusterType.lower() or 'both' in self.clusterType.lower() or 'combined' in self.clusterType.lower():
			self.allClusters = pd.concat([self.gc, self.oc])
			self.df = self.allClusters
			# List to iterate through for all cluster plots later
			self.cluster_dfs = [self.gc, self.oc]

	def getClusterStats(self, params_to_write, output_file = True):
		'''
		Method to produce statistics for different cluster params 
		6 in total:
			-Cluster Mass (Msol)
			-Half-Mass radius (pc)
			-Distance (kpc)
			-Metallicity
			-Velocity Dispersion (km/s)
			-Cluster Age (Myr)

		params_to_write: list-like of cluster parameters to generate statristics for
		output_file: (boolean) if True, cluster stats are written to .txt files in 'output' directory
		'''
		# Drawing param data from dataframes
		self.gcMass, self.ocMass = self.gc['mass[Msun]'], self.oc['mass[Msun]']  # Mass
		self.gcRhm, self.ocRhm = self.gc['rhm[pc]'], self.oc['rhm[pc]']  # Rhm
		self.gcDist, self.ocDist = self.gc['dist[pc]'], self.oc['dist[pc]']  # Distance
		self.gcAge, self.ocAge = self.gc['age[Myr]'], self.oc['age[Myr]']  # Age
		self.gcVD, self.ocVD = self.gc['sigma_v0_z[km/s]'], self.oc['sigma_v0_z[km/s]']  # Velocity Dispersion
		self.gcZ, self.ocZ = self.gc['[Fe/H]'], self.oc['[Fe/H]']  # Metallicity

		
		# Printing mean
		print('Cluster Stats: ' + '\n')
		for c in params_to_write:
			self.paramLabel = self.paramDict[c.lower()]
			self.unit = self.unitsDict[c.lower()]
			# Printing GC param valuues to screen (summed and averaged)
			print(f'Globular Clusters ({c}): ')
			print(f'The mean value for globular cluster {c} is: {np.mean(self.gc[self.paramLabel])}' + ' ' + self.unit)
			print(f'The total (summed across all clusters) value for globular cluster {c} is: {np.sum(self.gc[self.paramLabel])}' + ' ' + self.unit + '\n')
			# print('')
			# Printing OC param valuues to screen (summed and averaged)
			print(f'Open Clusters ({c})')
			print(f'The mean value for open cluster {c} is: {np.mean(self.oc[self.paramLabel])}' + ' ' + self.unit)
			print(f'The total (summed across all clusters) value for open cluster {c} is: {np.sum(self.oc[self.paramLabel])}' + ' ' + self.unit + '\n')

		# Also print out all cluster values

		# If desired, writing cluster stats to output txt files
		if output_file:
			print('Writing data to txt output file...')
			with open(os.path.join('output', f'{self.clusterType}-Stats.txt'), 'w') as f:
				for c in params_to_write:
					self.paramLabel = self.paramDict[c.lower()]
					self.unit = self.unitsDict[c.lower()]
					f.write(c.upper() + '\n')
					f.write('Globular Clusters: ' + '\n')
					f.write(f'The mean value for globular cluster {c} is: {np.mean(self.gc[self.paramLabel])}' + ' ' + self.unit  + '\n')
					f.write(f'The total (summed across all clusters) value for globular cluster {c} is: {np.sum(self.gc[self.paramLabel])}' + ' ' + self.unit + '\n')
					# f.write('')
					f.write('Open Clusters: ' + '\n')
					f.write(f'The mean value for open cluster {c} is: {np.mean(self.oc[self.paramLabel])}' + ' ' + self.unit + '\n')
					f.write(f'The total (summed across all clusters) value for open cluster {c} is: {np.sum(self.oc[self.paramLabel])}' + ' ' + self.unit + '\n')
					f.write('######################################## \n')
			# Writing file with both cluster types as well for each param
				# if combined:


	def clusterHist(self, params_to_plot, show_mean=False, save_fig=False):
		'''
		Method to generate histogram plots of different cluster params
		params_to_plot - list-like of paramters to geenrate histograms for

		'''
		# for d in cluster_dfs:
		for p in params_to_plot:
			# Param labels correspond to labels used in original .csv files
			# These are different than the keys passed to clusterHist (paramLabel includes units)
			paramLabel = self.paramDict[p.lower()]
			for d in self.cluster_dfs:  # if combined, will plot OCs and GCs in different colors
				# Fix legend labelins
				# if 'gc' in self.clusterType.lower() or 'globular' in self.clusterType.lower():
				plt.hist(d[paramLabel], bins=30, label=self.clusterType, alpha=0.5)

				plt.title(self.clusterType)
				plt.xlabel(paramLabel)
				plt.legend()
				# If user wants mean value showed on plot as a vertical line
				if show_mean:
					plt.vlines(np.mean(d[paramLabel]), 0, len(d), linestyle='--', label=f'mean:{clusterType}')
			# plt.show()

			if save_fig:
				path = os.path.join('plots', 'hists')
				plt.savefig(os.path.join(path, f'{self.clusterType}-{p}-hist.pdf'), format='pdf')

	def clusterScatter(self, params, save_fig=False):
		'''
		Method to generate 2D histogram between 2 cluster params
		Done instead of a scatter 
		params - list-like of cluster params (give param key strings)
		clusterType - used for plot title/labels

		'''
		# x and y data pulled from params key list
		print(self.df)
		x = self.df[self.paramDict[params[0]]]
		y = self.df[self.paramDict[params[1]]]

		plt.scatter(x,y)
		plt.xlabel(self.paramDict[params[0]])
		plt.ylabel(self.paramDict[params[1]])
		plt.title(self.clusterType)
		# plt.show()

		# If desired, saving figure to plots dir
		if save_fig:
			path = os.path.join('plots', 'scatter')
			plt.savefig(os.path.join(path, f'{self.clusterType}-{params[0]}-{params[1]}-scatter.pdf'), format='pdf')

	def clusterHist2d(self, params, save_fig=False):
		'''
		Method to generate 2D histogram between 2 cluster params
		Done instead of a scatter 
		params - list-like of cluster params (give param key strings)
		clusterType - used for plot title/labels

		'''
		# x and y data pulled from params key list
		x = self.df[self.paramDict[params[0]]]
		y = self.df[self.paramDict[params[1]]]

		plt.hist2d(x,y)
		plt.xlabel(self.paramDict[params[0]])
		plt.ylabel(self.paramDict[params[1]])
		plt.title(self.clusterType)
		# plt.show()

		# If desired, saving figure to plots dir
		if save_fig:
			path = os.path.join('plots', 'hists', 'hist2d')
			plt.savefig(os.path.join(path, f'{self.clusterType}-{params[0]}-{params[1]}-hist2D.pdf'), format='pdf')


# #############################################################################################################

# Test call of object/methods
# cs = clusterStats(['mass', 'age', 'rhm', 'dist', 'Z', 'sigma'], 'OC')
# cs.getClusterData()
# cs.getClusterStats(['mass', 'age', 'rhm', 'dist', 'Z', 'sigma'], True)
# # cs.clusterScatter(['dist', 'rhm'])
# # cs.clusterHist2d(['mass', 'age'], False)


# TODO: 
# - Add more plotting methods
# - Test output file writing

types = ['OC', 'GC', 'all']
params = ['mass', 'age', 'rhm', 'dist', 'z', 'sigma']

for t in types:
	print(f'Producing cluster data for {t}')

	C = clusterStats(params, t)
	C.getClusterData()
	C.clusterHist(params, False, True)
	for p in permutations(params,2):

		C.clusterScatter(p, True)
		C.clusterHist2d(p, True)

	










