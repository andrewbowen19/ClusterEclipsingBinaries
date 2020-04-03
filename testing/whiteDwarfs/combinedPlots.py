"""
Script to create combined histograms and scatter plots 
All 4 scenarios for each cluster type/subpopulation: 
i.e. recovered binaries for Globular Clusters: CN, CC, BC, BN all in one plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import permutations
# List and dict of labels and markerstyles for plots
markers = ['^', 'o', 's', '*']
labels = {'BN': 'Baseline No Crowding', 'BC':'Baseline Crowding', 'CN':'Colossus No Crowding', 'CC':'Colossus Crowding'}
alphas = np.arange(1.0,0.0,0.25)

def makePlottingArrays(path, files):
	'''
	Function to create combined scatterplot: diff colors/markers for all 4 scenarios
	Params:
	path - path to dubir and requisite files
	files - should be list of files in either Rec or Obs subdir (len(f) = 4)
	x and y labels - strings, pull from dataframe names
	title - can just pull from file name (e.g. rec-M10BN)
	scenario - pull from file name
	'''
	# Rewriting this function to do all the work of the below file walk
	m1 = []
	m2 = []
	r1 = []
	r2 = []
	scenarios = [] # observing scenario of that file (should be in order)
	for f in files:
		df = pd.read_csv(path + f)
		df.drop('Unnamed: 0', axis =1)

		m1.append(df['m1'])
		m2.append(df['m2'])
		r1.append(df['r1'])
		r2.append(df['r2'])

		# Generating appropriate observing scenario string
		if 'G' in f or 'O' in f:
			scenario = f[5:7]
			scenarios.append(scenario)
		elif 'M10' in f or 'M67' in f:
			scenario = f[7:9]
			scenarios.append(scenario)


		

	return {'m1': m1, 'm2': m2, 'r1': r1, 'r2':r2}, scenarios


def multiScatter(data_x, data_y, xlabel, ylabel, title,scenario_list):
	'''
	Method to actually make the scatter plots
	data_x and data_y - nested lists of len 4 (1 entry for each scenario in a dir)
	scenario_list - ordered list of observing scenarios

	'''
	f,ax = plt.subplots()

	for i in range(0,len(data_x)):
		ax.scatter(data_x[i], data_y[i], marker = markers[i], alpha = 0.5, label = scenario_list[i])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	ax.legend()
	# plt.show()
	f.savefig(f'./plots/combinedPlots/scatter/{title}-{xlabel}{ylabel}-scatter.pdf')
	print('Scatter plots made!')

# #################################################################################################################

# Looping through WD binary file tree to make these combined scatter plots and hists
for root, dirs, files in os.walk('./wd_output/', topdown = True):
	for d in dirs:
		print('working dir: ', d)
		# haveObs = False
		# haveRec = False
		# m1Obs = []
		# m2Obs = []
		# r1Obs = []
		# r2Obs = []

		# m1Rec = []
		# m2Rec = []
		# r1Rec = []
		# r2Rec = []

		nameslices = []
		scenarioSlices = []

		# Getting correct slices for different cluster paths
		if 'GlobularClusters' in root or 'OpenClusters' in root:
			nameslices = [0,7]
			scenarioSlices = [5,7]
		if 'm10' in root or 'm67' in root:
			sameslices = [0,9]
			scenarioSlices = [7,9]

		# Only want 2 subdirs
		if 'Obs' in d or 'Rec' in d:

			print(d)
			f = os.listdir(root + '/' + d)
			print(f)

			# Picking out correct name slices for plot tiles (i.e. 'rec-G' for all recovered binaries in GCs)
			for name in f:
				if 'G' in name or 'O' in name:
					plot_title = name[0:5]
				if 'M10' in name or 'M67' in name:
					plot_title = name[0:7]

			# Reading in data files and creating nested lists for plotting
			dat, obs_scen = makePlottingArrays(root + '/' + d + '/',f) # makes correctly formatted nested lists for each Obs and Rec subdir

			# Looping through all param permutations (between m1,m2,r1,r2)
			for p in permutations(dat, 2):
				# x and ya labels
				x_param = p[0]
				y_param = p[1]
				
				# Making scatter plots for all combinations of mass-radius params
				multiScatter(dat[x_param], dat[y_param], p[0], p[1], plot_title, obs_scen)




