"""
Script to create combined histograms and scatter plots
All 4 scenarios for each cluster type/subpopulation:
i.e. recovered or observed binaries for Globular Clusters: CN, CC, BC, BN all in one plot
Using this to generate scatter plots for WD binary pairs
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import permutations

# List and dict of labels and markerstyles for plots
markers = ['^', 'o', 's', '*']
labels = {'BN': 'Baseline No Crowding', 'BC':'Baseline Crowding', 
			'CN':'Colossus No Crowding', 'CC':'Colossus Crowding'}

colors = ['#a6cee3' , '#1f78b4', '#b2df8a' , '#33a02c']

def makePlottingArrays(path, files):
	'''
	Function to create combined scatterplot: diff colors/markers for all 4 scenarios
	Params:
	path - path to dubir and requisite files
	files - should be list of files in either Rec or Obs subdir (len(f) = 4)
	'''

	# Rewriting this function to do all the work of the below file walk
	m1 = []
	m2 = []
	r1 = []
	r2 = []
	scenarios = [] # observing scenario of that file (should be in order)
	p = []
	for f in files:
		df = pd.read_csv(path + f)
		df.drop('Unnamed: 0', axis =1)

		# Appending data to df
		m1.append(df['m1'])
		m2.append(df['m2'])
		r1.append(df['r1'])
		r2.append(df['r2'])

		p.append(df['p'])

		# Generating appropriate observing scenario string
		if 'G' in f or 'O' in f:
			scenario = f[5:7]
			scenarios.append(scenario)
		elif 'M10' in f or 'M67' in f:
			scenario = f[7:9]
			scenarios.append(scenario)


		

	return {'m1': m1, 'm2': m2, 'r1': r1, 'r2':r2}, scenarios, p

def multiScatter(data_x, data_y, xlabel, ylabel, title,scenario_list):
	'''
	Method to actually make the scatter plots
	data_x and data_y - nested lists of len 4 (1 entry for each scenario in a dir)
	scenario_list - ordered list of observing scenarios

	'''
	f,ax = plt.subplots()
	# Plotting for all 4 viewing scenarios
	for i in range(0,len(data_x)):
		ax.scatter(data_x[i], data_y[i], marker = markers[i], alpha = 0.5, c = colors[i],label = scenario_list[i])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	ax.legend()

	f.savefig(f'./plots/combinedPlots/scatter/{title}-{xlabel}{ylabel}-scatter.pdf')
	print('Scatter plots made!')

def multiHist(data, title, xlabel, xmin, xmax, scenario_list):
	'''
	Function to create layered histograms among all 4 params
	Will likely use for period histograms among all 4 scenarios
	Same arguments as multiScatter
	'''
	f,ax=plt.subplots()
	for i in range(0, len(data)):
		ax.hist(data[i], range=[xmin, xmax],  histtype='step', color=colors[i], label=scenario_list[i])
	ax.set_xlabel(xlabel)
	ax.set_title(title)
	ax.legend()
	f.savefig(f'./plots/combinedPlots/hists/{title}-{xlabel}Hist.pdf')
	print(f'Combined {xlabel} histgram made!')

# #################################################################################################################

# Looping through WD binary file tree making combined scatter plots and hists
for root, dirs, files in os.walk('./wd_output/', topdown = True):
	for d in dirs:
		print('working dir: ', d)

		nameslices = []
		scenarioSlices = []

		# Getting correct slices for different cluster paths
		if 'GlobularClusters' in root or 'OpenClusters' in root:
			nameslices = [0,7]
			scenarioSlices = [5,7]
		if 'm10' in root or 'm67' in root:
			sameslices = [0,9]
			scenarioSlices = [7,9]

		# Only want 2 subdirs - observed and recovered
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
			dat, obs_scen, periods = makePlottingArrays(root + '/' + d + '/',f) # makes correctly formatted nested lists for each Obs and Rec subdir

			# Making combined period histograms
			multiHist(periods, plot_title, 'p', 0., 50., obs_scen)

			# Making mass, radius combined hists
			multiHist(dat['m1'], plot_title, 'm1', 0., 0.6, obs_scen) # m1
			multiHist(dat['m2'], plot_title, 'm2', 0., 0.6, obs_scen) # m2
			multiHist(dat['r1'], plot_title, 'r1', 0., 0.8, obs_scen) # r1
			multiHist(dat['r2'], plot_title, 'r2', 0., 0.8, obs_scen) # r2

			# Looping through all param permutations (between m1,m2,r1,r2) -- making scatter plots
			for p in permutations(dat, 2):
				# x and ya labels
				x_param = p[0]
				y_param = p[1]
				
				# Making scatter plots for all combinations of mass-radius params
				multiScatter(dat[x_param], dat[y_param], p[0], p[1], plot_title, obs_scen)




