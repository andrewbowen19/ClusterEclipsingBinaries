'''
Corner plotting script to produce corner plots with multiple colors for different cluster types
Still loops through current file path to create corner plots for each observing scenario
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
def corner_plot(DataFrame1, DataFrame2, popType, name, contours = False):
	'''
	Function to make a corner plot with our different binary parameters -- using corner.py
	Creates corner plot for two different cluster sub-pops (OC vs GC)
	df1 - should be dataframe with values (all binary params) (each columnn is one element) (OCs)
	df2 - same structure as df1, but for GCs
	popType - binary subpopulation type (All, Obs, Rec) 
	name - name of scenario (Open/Globular/Long; Baseline/Colossus; Crowding/NoCrowding)
	contours - set as True if contours overlaid on plots desired, default: False
	'''
	# Plot with Contours
	if contours == True:
		print('Making corner plots with contours...')
		df1 = DataFrame1
		df2 = DataFrame2
		# Plotting first dataframe - OCs
		f1 = corner.corner(DataFrame1, color = 'r', labels = df1.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, scale_hist = True, hist2d_kwargs = {'plot_density':True})  # , scale_hist = True)
		# Plotting secodn dataframe - GCs
		corner.corner(DataFrame2, color = 'b', labels = df2.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, fig = f1, hist2d_kwargs = {'plot_density':True})
		f1.suptitle(popType + '-' + name, fontsize=24)
		plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/', f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner contour plots made!')

	# No contours
	elif contours == False:
		print('Making corner plots...')
		df1 = DataFrame1
		df2 = DataFrame2
		f1 = corner.corner(df1, color = 'r', labels = df1.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = False, title_kwargs={"fontsize": 28})# , scale_hist = True)
		corner.corner(df2, color = 'b', labels = df2.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = False, title_kwargs={"fontsize": 28}, scale_hist = True, fig = f1)
		f1.suptitle(popType + '-' + name, fontsize=24)
		plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/', f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner plots made!')

	print('On to the next!')

	
# ########################################################################################################

ocPath = os.path.join('clusters', 'OpenClusters', 'baseCrowd', 'data')
gcPath = os.path.join('clusters', 'GlobularClusters', 'baseCrowd', 'data')

for ocFile, gcFile in zip(os.listdir(ocPath), os.listdir(gcPath)):
	
	# Ignoring wd file dsubdirs
	if '.csv' in ocFile and '.csv' in gcFile:
		print('File name:', ocFile, gcFile)

		# Reading in binary data for all 3 sub-pops
		ocDat = pd.read_csv(os.path.join(ocPath, ocFile), header = 0)
		gcDat = pd.read_csv(os.path.join(gcPath, gcFile), header = 0)

		# Changing columns/period to log
		ocDat['p'] = np.log10(ocDat['p']) # converting period values from days to log-days
		gcDat['p'] = np.log10(gcDat['p']) # converting period values from days to log-days
		ocDat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$',
										'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r']
		gcDat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$',
										'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r'] 
		# print('Binary Pop data for OCs, GCs: ', ocDat, '\n', gcDat)

		# Making multi corner plots - OCs in red, GCs in blue
		corner_plot(ocDat, gcDat, ocFile[0:3], ocFile[4:7] + '-' + gcFile[4:7], True)


# Looping through correct files in our trees and making plots
# for root, dirs, files in os.walk('.', topdown = True):

# 	for name in files:

# 		print('NAME: ',name)
# 		if '-histData.csv' in name:

# 			# Reading in an renaming dataframe columns for plotting - doing once and outside of our function call
# 			dat = pd.read_csv(os.path.join(root, name), header = 0)
# 			dat['p'] = np.log10(dat['p']) # converting period values from days to log-days
# 			dat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$', 'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r']


# 			print('Data read in...')
			
# 			# Making corner plots for every scenario -- add to other makeHists scripts (M67, OCs and GCs)
# 			if 'GlobularClusters' in root or 'OpenClusters' in root:
# 				corner_plot(dat, name[0:3], name[4:7], True) # WITH contours -- Change name slice to [4:7] for OCs and GCs (less characters)
# 				corner_plot(dat, name[0:3], name[4:7], False) # NO contours

# 			elif 'm10' in root or 'm67' in root:
# 				corner_plot(dat, name[0:3], name[4:9], True) # WITH contours -- Change name slice to [4:7] for OCs and GCs (less characters)
# 				corner_plot(dat, name[0:3], name[4:9], False) # NO contours

# print('All done!')





