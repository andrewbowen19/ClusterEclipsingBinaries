"""
Script to create combined histograms and scatter plots 
All 4 scenarios for each cluster type/subpopulation: 
i.e. recovered binaries for Globular Clusters: CN, CC, BC, BN all in one plot
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# List and dict of labels and markerstyles for plots
markers = ['^', 'o', 's', '*']
labels = {'BN': 'Baseline No Crowding', 'BC':'Baseline Crowding', 'CN':'Colossus No Crowding', 'CC':'Colossus Crowding'}
alphas = np.arange(1.0,0.0,0.25)

def makeCombinedScatterPlot(data_x, data_y, xlabel, ylabel, title, scenario):
	'''
	Function to create combined scatterplot: diff colors/markers for all 4 scenarios
	Params:
	data - listlike object containing all 4 scenarios (BC, BN, CC, CN) for each cluster type and recovery pop (GC/OC/m10/m67; rec/obs)
		Should be a nested list of 4 scenario lists of the data we want to plot
	x and y labels - strings, pull from dataframe names
	title - can just pull from file name (e.g. rec-M10BN)
	scenario - pull from file name
	'''
	# Selecting scenario labels based on plot title
	if 'G' in title or 'O' in title:
		label = title[]

	f,ax = plt.subplots()

	for i in range(0,len(data_x)-1):
		ax.scatter(data_x[i], data_y[i], marker = markers[i], alpha = 0.5, label = labels[scenario])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	ax.legend()
	plt.show()
	f.savefig(f'./plots/combinedPlots/scatter/{title}-{xlabel}{ylabel}-scatter.pdf')

# Test call for function
m10rec = ['rec-M10BC-wdBinaries.csv', 'rec-M10BN-wdBinaries.csv', 'rec-M10CC-wdBinaries.csv', 'rec-M10CN-wdBinaries.csv']

m1 = []
m2 = []
for file in m10rec:
	dat = pd.read_csv('./wd_output/m10/RecM10/'+file)
	dat = dat.drop('Unnamed: 0', axis =1)
	m1.append(dat['m1'])
	m2.append(dat['m2'])

print(len(m1), range(len(m1)))
xx = makeCombinedScatterPlot(m1, m2, 'm1', 'm2', 'rec-M10', labels[])

# Looping through WD binary file tree to make these combined scatter plots and hists
for root, dirs, files in os.walk('./wd_output/', topdown = True):


	for d in dirs:
		print('foo', d)
		haveObs = False
		haveRec = False
		m1Obs = []
		m2Obs = []
		r1Obs = []
		r2Obs = []

		m1Rec = []
		m2Rec = []
		r1Rec = []
		r2Rec = []

		# Only want 2 subdirs
		if 'Obs' in d or 'Rec' in d:

			for name in files:
				print(root, d, name)
				if 'wdBinaries.csv' in name:
					# print('Root, Name: ', root, name)

					if 'obs' in name:
						haveObs = True
						datObs = pd.read_csv(root + '/' + name)
						datObs = datObs.drop('Unnamed: 0', axis =1)
						# Appending data to lists
						m1Obs.append(datObs['m1'])
						m2Obs.append(datObs['m2'])
						r1Obs.append(datObs['r1'])
						r2Obs.append(datObs['r2'])


					elif 'rec' in name:
						haveRec = True
						datRec = pd.read_csv(root + '/' + name)
						datRec = datRec.drop('Unnamed: 0', axis =1)
						# Appending data to lists
						m1Rec.append(datRec['m1'])
						m2Rec.append(datRec['m2'])
						r1Rec.append(datRec['r1'])
						r2Rec.append(datRec['r2'])

					# print(m1Rec, m1Obs)

					if haveObs:
						print('making Obs scatter plots...')
						makeCombinedScatterPlot(m1Obs, m2Obs, 'm1', 'm2',name[0:7], name[5:8])
						
					if haveRec:
						print('making Rec scatter plot...')
						makeCombinedScatterPlot(m1Rec, m2Rec, 'm1', 'm2',name[0:7], name[5:8])








