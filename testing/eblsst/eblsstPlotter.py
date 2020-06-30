# Script to create plots from output eblsst files
# These files are grabbed from different observing scenarios
# They are produced from the analyseClusterLISA script located in the flex directory
# Let's plot!

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


class eblsstPlotter(object):
	'''
	Class to create histograms for 
	'''
	def __init__(self, DataFrame, clusterType, strategy, crowding, param):
		self.DataFrame = DataFrame
		self.clusterType = clusterType
		self.strategy = strategy
		self.param = param  # binary parameter to be plotted
		# Filter plotting colors '<filter name>': '<matplotlib color>'
		# More info here: https://community.lsst.org/t/lsst-filter-profiles/1463
		self.filterColors = {'u': '#8B00FF', 'g': '#0000FF',
							 'r': '#00FF00', 'i': '#FFFF00',
							 'z': '#FF7F00', 'y': '#FF0000',
							 'ugrizy': '#000000'}
		# Dictionary including parameter keys (from file names to more descriptive labels)
		self.paramDict = {'lp': 'log-period', 'd': 'dist',
						  'm1': 'm1', 'q': 'q',
						  'mag': 'mag', 'e': 'e', 'r': 'r'}
		self.crowding = crowding

	def makeIndHist(self, data, filter, show_mean=False, save_fig=True):
		'''
		Method to generate histograms of binaries for individual Vera Rubin filters
		data: list-like containing binaries in each bin (should be a 'hist' column from eblsst files)
		bins: bins for hist, use 'binsEdges' columnn from csv files
		'''
		print(f'Making Individual histogram for {self.clusterType}-{self.strategy}-{self.param}...')
		print(f'filter: {filter}')
		f, ax = plt.subplots(figsize=(8, 5))
		ax.hist(data, bins=self.DataFrame['binEdges'], histtype='step', color=self.filterColors[filter])
		ax.set_title(self.clusterType + '-' + self.strategy + self.crowding + '; ' + filter)
		ax.set_xlabel(self.param)#.replace('EBLSST_', ''))

		if show_mean:
			ax.vline(np.mean(data), linestyle='--', label='Distribution mean')

		if save_fig:
			fileName = f'{self.clusterType}-{self.strategy}-{self.param}-{filter}-histInd.pdf'
			path_to_fig = os.path.join('.', 'plots', self.clusterType, self.strategy + self.crowding, 'ind', fileName)
			f.savefig(path_to_fig, format='pdf')
		# f.close()
		print('Individual histogram complete!')



	def makeCombinedHist(self, filters, save_fig=True):
		'''
		Method to make histogram of binaries including filters
		data: just pass dataframe with filter columns to this
		filters: list of strings containinf filter keys
		'''
		multiData = self.DataFrame.drop(['binEdges', 'histAll', 'histObs', 'allhistRec'], axis=1)
		
		fig, ax = plt.subplots(figsize=(8, 5))

		for f in filters:
			ax.hist(multiData[f+'_histRec'], bins=self.DataFrame['binEdges'],
							  histtype='step', linestyle='--', color=self.filterColors[f], label=f, alpha=0.7)
		ax.legend()
		ax.set_title(ax.set_title(self.clusterType + '-' + self.strategy + self.crowding + '; recovered'))
		ax.set_xlabel(self.paramDict[self.param])  # replacing parameter label from file name w param key (in dict above)
		ax.set_ylabel(r'$N_{bin}$')

		if save_fig:
			filterStr = ''.join(filters)
			fileName = f'{self.clusterType}-{self.strategy}-{self.param}-{filterStr}-histCombined.pdf'
			path_to_fig = os.path.join('.', 'plots', self.clusterType, self.strategy + self.crowding, 'combined', fileName)
			fig.savefig(path_to_fig)


# #########################################################################################################################

# Looping through all observing scenarios
# path_to_clusters = os.path.join('.', 'clusters')
# for root, dirs, files in os.walk(path_to_clusters, topdown=True):
# 	print(root)
# 	for 





















# ########################################################################################################################

# Old path loop: only did GBC scenario
# test path to play around with plotting, will loop through file tree later
# test_path = os.path.join('clusters', 'GlobularClusters', 'baseCrowd', 'eblsst_files')
# save_path = test_path.replace('clusters', 'plots')


# filtersLSST = ['u', 'g', 'r', 'i', 'z', 'y']

# # looping through one dir for now
# for file  in os.listdir(test_path):
#     print(file)
#     df = pd.read_csv(os.path.join(test_path, file), header=0)
#     print(df)
#     binParam = file.replace('EBLSST_', '').replace('hist.csv', '')
#     eb = eblsstPlotter(df, 'GlobularClusters', 'base', 'Crowd', binParam)
#     eb.makeCombinedHist(filtersLSST, True)

#     for c in df.columns:

#     	# Getting filter key from column titles (if an individual filter col)
#     	# If there is an individual filter column, it's recovered
#     	if '_histRec' in c:
#     		filt = c.replace('_histRec', '')
#     		eb.makeIndHist(df[c],  filt, False, True)
#     	else:
#     		eb.makeIndHist(df[c], 'ugrizy', False, True)





# binEdges,histAll,histObs,u_histRec,g_histRec,r_histRec,i_histRec,z_histRec,y_histRec,allhistRec

