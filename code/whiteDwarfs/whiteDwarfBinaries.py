"""
Script for finding WD-WD recovered binary pairs
pulling data from our WD-histData.csv files
Mostly care about rec and obs subpopulations
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Maybe write this as a class for practice
class whiteDwarfBinaries(object):

	def __init__(self, path, filename):
		self.path = './'
		self.filename = '*WD-histData.csv'

	def wdMRrelation(self, mass):
		'''
		Function for mass-radius relationship for MS stars, 
		WDs will have radii smaller than the predicted R from this relation
		Source: http://personal.psu.edu/rbc3/A534/lec18.pdf
		'''
		xi = 0
		for m in mass:
			if (m < 1.0):
				xi = 0.8
				radius = m ** xi
				# print('radius: ', radius)
				return radius

			elif m >= 1.0:
				xi = 0.57
				radius = m ** xi
				# will not have negative radii values
				return 0.0

	def fileReadIn(self, path, filename):
		'''
		Function to read in WD csv files
		Should all reside in same subdirectory structure: <viewing scenario>/data/wd/
		'''
		self.path = path
		self.filename = filename
		if 'WD-histData.csv' in self.filename:
			self.df = pd.read_csv(self.path + self.filename) #, names = ['p','m1','m2','r1','r2','e','i','appMagMean_r'])
			self.df = self.df.drop('Unnamed: 0', axis =1) # dropping index column from csv write

			return self.df

	def scatterPlot(self, x,y, title, xlabel, ylabel, saveFig = False):
		'''
		Function to make a scatter plot between two variables, (x and y) 
		Provide labels and titles as needed
		'''
		f,ax = plt.subplots()
		ax.scatter(x,y)
		ax.set_xlabel(xlabel, fontsize = 16)
		ax.set_ylabel(ylabel, fontsize = 16)
		ax.set_title(title, fontsize = 18)

		# plt.show()

		if saveFig == True:
			# labels to include in saved figure name
			xparam = xlabel[1:4].replace('_', '')
			yparam = ylabel[1:4].replace('_', '')

			# Saving to correct cluster type subdir in /plots/whiteDwarf/ pairs dir
			# Globular Cluster
			if 'G' in title:
				f.savefig(f'./plots/whiteDwarfPairs/GlobularClusters/{title}-WDscatter-{xparam}-{yparam}.pdf')

			# Open Clusters
			elif 'O' in title:
				f.savefig(f'./plots/whiteDwarfPairs/OpenClusters/{title}-WDscatter-{xparam}-{yparam}.pdf')

			# M10
			elif 'M10' in title:
				f.savefig(f'./plots/whiteDwarfPairs/m10/{title}-WDscatter-{xparam}-{yparam}.pdf')

			# M67
			elif 'M67' in title:
				f.savefig(f'./plots/whiteDwarfPairs/m67/{title}-WDscatter-{xparam}-{yparam}.pdf')

	def pHist(self, DataFrame, scenario, saveFig = True):
		'''
		Method to make period histogram (periods in units of days)
		'''
		f,ax = plt.subplots()
		ax.hist(DataFrame['p'], bins = 10)
		ax.set_xlabel('period (days)', fontsize = 16)
		ax.set_title(scenario, fontsize = 18)

		# Saving figures to plots/whiteDwarfPairs subdir
		if saveFig == True:
			# Globular Cluster
			if 'G' in scenario:
				f.savefig(f'./plots/whiteDwarfPairs/GlobularClusters/{scenario}-WDhist-p.pdf')

			# Open Clusters
			elif 'O' in scenario:
				f.savefig(f'./plots/whiteDwarfPairs/OpenClusters/{scenario}-WDhist-p.pdf')

			# M10
			elif 'M10' in scenario:
				f.savefig(f'./plots/whiteDwarfPairs/m10/{scenario}-WDhist-p.pdf')

			# M67
			elif 'M67' in scenario:
				f.savefig(f'./plots/whiteDwarfPairs/m67/{scenario}-WDhist-p.pdf')


	def findWhiteDwarfPairs(self, DataFrame, scenario, plots = True, save_csv = True):
		'''
		Function to select white dwarf binary pairs (both binaries meet WD criteria we used in analyse plots)
		Params: DataFrame - dataframe object to be passed to our function
		plots - whether or not user wants to generate plots
		scenario - which combo of viewing factors (ex. rec-GBC or all-M67CN)
		'''

		# Only selecting binary pairs that both meet selection criteria
		self.wds = DataFrame.loc[((DataFrame['m1'] < 0.6) & (DataFrame['m2'] < 0.6))  
			& ((DataFrame['r1'] < self.wdMRrelation(DataFrame['m1'])) & (DataFrame['r2'] < self.wdMRrelation(DataFrame['m2'])))]
		
		if plots == True:
			# histogram of periods (no need for log)
			self.pHist(self.wds, scenario, True)

			# Scatter plot of primary masses and radii
			self.scatterPlot(self.wds['m1'], self.wds['r1'],scenario, '$m_1$ $(M_{\odot})$','$r_1$ $(R_{\odot})$', True)

			# Scatter plot of secondary masses and radii
			self.scatterPlot(self.wds['m2'], self.wds['r2'], scenario, '$m_2$ $(M_{\odot})$','$r_2$ $(R_{\odot})$', True)

			# Mass scatters
			self.scatterPlot(self.wds['m1'], self.wds['m2'], scenario, '$m_1$ $(M_{\odot})$','$m_2$ $(M_{\odot})$', True)

			# Radius scatter plots
			self.scatterPlot(self.wds['r1'], self.wds['r2'], scenario, '$r_1$ $(R_{\odot})$','$r_2$ $(M_{\odot})$', True)

		if save_csv == True:
			self.wds.to_csv(f'./wd_output/{scenario}-wdBinaryPairs.csv')

		# return wds


# Test call of readin function

# wdb = whiteDwarfBinaries('./m10/baseNoCrowdM10/data/wd/', 'rec-M10BN-WD-histData.csv')
# dat = wdb.fileReadIn('./m10/baseNoCrowdM10/data/wd/', 'rec-M10BN-WD-histData.csv')
# plots = wdb.findWhiteDwarfPairs(dat, 'rec-M10BN', True, True)


# Periods in our WD data files are NOT in the log

# Now we'll loop through our file tree to create these plots for 
for root, dirs, files in os.walk('./clusters/', topdown = True):
	for name in files:

		if ('WD-histData.csv' in name) and ('rec' in name or 'obs' in name):
			print('ROOT, NAME : ', root, name)
			print('Found file, reading in data...')
			wdb = whiteDwarfBinaries(root + '/', name)

			data = wdb.fileReadIn(root + '/', name)
			print(data)
			plotMe = wdb.findWhiteDwarfPairs(data, name.replace('-WD-histData.csv', ''), True, True)
			print('Plots made!')
			print('')











