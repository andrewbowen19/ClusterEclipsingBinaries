# Script to create plots from output eblsst files
# These files are grabbed from different observing scenarios
# They are produced from the analyseClusterLISA script located in the flex directory
# Let's plot!

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# test path to play around with plotting, will loop through file tree later
test_path = os.path.join('clusters', 'GlobularClusters', 'baseCrowd', 'eblsst_files')

class eblsstPlotter(object):
	'''
	Class to create histograms for 
	'''
	def __init__(self, DataFrame, clusterType, strategy, param):
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

	def makeIndHist(self, data, bins, filter, show_mean=False, save_fig=True):
		'''
		Method to generate histograms of binaries for individual Vera Rubin filters
		data: list-like containing binaries in each bin (should be a 'hist' column from eblsst files)
		bins: bins for hist, use 'binsEdges' columnn from csv files
		'''
		print(f'Making Individual histogram for {self.clusterType}-{self.strategy}-{self.param}...')
		print(f'filter: {filter}')
		f, ax = plt.subplots(figsize=(8, 5))
		ax.hist(data, bins=bins, histtype='step', color=self.filterColors[filter])
		ax.set_title(self.clusterType + '-' + self.strategy + '; ' + filter)
		ax.set_xlabel(self.param.replace('EBLSST_', ''))

		if show_mean:
			ax.vline(np.mean(data), linestyle='--', label='Distribution mean')

		if save_fig:
			fileName = f'{self.clusterType}-{self.strategy}-{self.param}-{filter}-histInd.pdf'
			path_to_fig = os.path.join('.', 'plots', fileName)
			f.savefig(path_to_fig, format='pdf')
		# f.close()
		print('Individual histogram complete!')



	# def makeCombinedHist(self, data, filters):
	# 	'''
	# 	Method to make histogram of binaries across filters
	# 	'''

# looping through one dir for now
for f in os.listdir(test_path):
    print(f)
    df = pd.read_csv(os.path.join(test_path, f), header=0)
    print(df)
    eb = eblsstPlotter(df, 'GC', 'base', f[0:9])
    for c in df.columns:

    	# Getting filter key from column titles (if an individual filter col)
    	# If there is an individual filter column, it's recovered
    	if '_histRec' in c:
    		filt = c.replace('_histRec', '')
    		eb.makeIndHist(df[c], df['binEdges'], filt, False, True)
    	else:
    		eb.makeIndHist(df[c], df['binEdges'], 'ugrizy', False, True)




