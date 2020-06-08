'''
Corner plotting script -- on Quest
Flex repo version, loops through current file path to create corner plots for each observing scenario
Updated for flex file paths (even for all observing scenarios)
Will create corner plots for all binary sample (i.e. not just WDs)
Corner plotting documentation: https://corner.readthedocs.io/en/latest/
'''

import corner
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg') # uncomment if running on Quest
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os

# Maybe put this in its own script for runtime purposes
def corner_plot(DataFrame, popType, name, contours = False):
	'''
	Function to make a corner plot with our different binary parameters -- using corner.py
	df - should be dataframe with values (all binary params) (each columnn is one element)
	popType - binary subpopulation type (All, Obs, Rec)
	name - name of scenario (Open/Globular/Long; Baseline/Colossus; Crowding/NoCrowding)
	contours - set as True if contours overlaid on plots desired, default: False
	'''
	# Plot with Contours
	if contours == True:
		print('Making corner plots with contours...')
		df = DataFrame
		f = corner.corner(df, labels = df.columns, label_kwargs={"fontsize":18}, bins = 20, plot_contours = True, title_kwargs={"fontsize": 28})
		f.suptitle(popType + '-' + name, fontsize=24)
		f.show()
		f.savefig(os.path.join('./plots/corner_plots/contours/', f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner contour plots made!')

	# No contours
	elif contours == False:
		print('Making corner plots...')
		df = DataFrame
		f = corner.corner(df, labels = df.columns, label_kwargs={"fontsize":18}, bins = 20, plot_contours = False, title_kwargs={"fontsize": 28})
		f.suptitle(popType + '-' + name, fontsize = 24)
		f.show()
		f.savefig(os.path.join('./plots/corner_plots/', f'{popType}-{name}-cornerPlot.pdf'))
		plt.close()
		print('Corner plots made!')

	print('On to the next!')

	
# ########################################################################################################

# Looping through correct files in our trees and making plots
for root, dirs, files in os.walk('.', topdown = True):

	for name in files:

		print('NAME: ',name)
		if '-histData.csv' in name:

			# Reading in an renaming dataframe columns for plotting - doing once and outside of our function call
			dat = pd.read_csv(os.path.join(root, name), header = 0)
			dat['p'] = np.log10(dat['p']) # converting period values from days to log-days
			dat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$', 'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r']


			print('Data read in...')
			
			# Making corner plots for every scenario -- add to other makeHists scripts (M67, OCs and GCs)
			if 'GlobularClusters' in root or 'OpenClusters' in root:
				corner_plot(dat, name[0:3], name[4:7], True) # WITH contours -- Change name slice to [4:7] for OCs and GCs (less characters)
				corner_plot(dat, name[0:3], name[4:7], False) # NO contours

			elif 'm10' in root or 'm67' in root:
				corner_plot(dat, name[0:3], name[4:9], True) # WITH contours -- Change name slice to [4:7] for OCs and GCs (less characters)
				corner_plot(dat, name[0:3], name[4:9], False) # NO contours

print('All done!')





