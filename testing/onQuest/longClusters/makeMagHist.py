# Script to make magnitude hists for long clusters (maybe grep for Glob/Open)

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


def makeMagHist(path, filename):


	df = pd.read_csv(path + filename, skiprows = 2) # skipping header rows (cluster data)
	clusterName = filename.replace('_output_file.csv', '')
	cluster = None
	strat = None
	crowd = None

	# Defining scenario
	if 'm67' in path:
		cluster = 'M67'
	elif 'm10' in path:
		cluster = 'M10'

	if 'colossus' in path:
		strat = 'colossus'
	elif 'colossus' not in path:
		strat = 'baseline'

	if 'withCrowding' in path:
		crowd = 'Crowd'
	if 'withCrowding' not in path:
		crowd = 'noCrowd'

	f,ax = plt.subplots()
	ax.hist(df['appMagMean_r'].loc[np.where(df['appMagMean_r'] != -999)]) # leaving out bad mags (a lot of them)
	ax.set_title(clusterName + ' ' + strat + ' ' + crowd)
	ax.set_xlabel('appMagMean_r') 
	plotname =  clusterName + strat + '-' +  crowd + '-magHist.pdf'
	f.savefig('./magPlots/' + cluster + '/' + plotname)
	# plt.show()


# making mag plots for long clusters
for root, dirs, files in os.walk('.', topdown = True):
	for name in files:

		if '_output_file.csv' in name:
			print(root, name)
			makeMagHist(root + '/', name)














