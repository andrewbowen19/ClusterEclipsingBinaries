'''
Script to run the eblsstPloter (located in this directory) for each observing scenario
The scenarios include every permutation of the following:
Open/Globular
baseline/colossus
crowding/no-crowding
Utilizes methods from eblsstPlotter to generate histograms for each scenario/binary param
*Note: in our nomenclature, OS refers to the 'observing scenario'. Sorry for any confusion.
'''

from eblsstPlotter import eblsstPlotter
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re


# dictionary including path folder names and the respective OpSim strategy and crowding scenario
scenarioDict = {'baseCrowd': ['base', 'Crowd'], 'baseNoCrowd': ['base', 'NoCrowd'],
				'colCrowd': ['col', 'Crowd'], 'colNoCrowd': ['col', 'NoCrowd']}
filtersLSST = ['u', 'g', 'r', 'i', 'z', 'y']

# Looping through all observing scenarios
path_to_clusters = os.path.join('.', 'clusters')
for root, dirs, files in os.walk(path_to_clusters, topdown=True):
	# print('ROOT: ', root)
	for d in dirs:
		# Grabbing data files from eblsst_files sub dirs (1 for each OS)
		if 'eblsst_files' in d:
			scenario = root.replace('./clusters/', '')
			print('Here\'s the sitch:', scenario)
			path_to_files = os.path.join(root , d)

			# Defining observing scenario (OC/GC, base/col, crowd/NoCrowd)
			cluster = scenario.split('/')[0]
			stratCrowd = scenario.split('/')[1]
			strategyOpSim = scenarioDict[stratCrowd][0]  # either baseline or colossus
			crowding = scenarioDict[stratCrowd][1]
			print(f'Making plots for following OS: {cluster}, {strategyOpSim}, {crowding}...')
			
			# Generating plots for each binary parameter file in working OS
			for file in os.listdir(path_to_files):

				df = pd.read_csv(os.path.join(path_to_files, file), header=0)
				binParam = file.replace('EBLSST_', '').replace('hist.csv', '')
				eb = eblsstPlotter(df, cluster, strategyOpSim, crowding, binParam)
				eb.makeCombinedHist(filtersLSST, True)
				

				for c in df.columns:

					# Getting filter key from column titles (if an individual filter col)
					# If there is an individual filter column, it's recovered
					if '_histRec' in c:
						filt = c.replace('_histRec', '')
						eb.makeIndHist(df[c],  filt, False, True)
					else:
						eb.makeIndHist(df[c], 'ugrizy', False, True)


print("All done!")








