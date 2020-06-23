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
matplotlib.use('Agg') # uncomment if running on Quest
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

		# Setting up better range values for either OCs or GCs
		ranges = [(-0.5, 2.5), (0., 15.), (0, 15.), (0., 8.), (0., 8.), (0.,1.), (45.,135.), (15.,24.)]

		# Plotting first dataframe - OCs
		f1 = corner.corner(DataFrame1, color = 'r', labels = df1.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, range = ranges)#, scale_hist = True, hist2d_kwargs = {'plot_density':True})  # , scale_hist = True)
		# Plotting secodn dataframe - GCs
		corner.corner(DataFrame2, color = 'b', labels = df2.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, fig = f1, range = ranges)# hist2d_kwargs = {'plot_density':True})
		f1.suptitle(popType + '-' + name, fontsize=24)
		#plt.show()
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
		#plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/', f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner plots made!')

	print('On to the next!')


def makeMultiCornerPlot(ocPath):
	'''
	Function to generate multi-colored corner plot
	Creates plot for an individual observing scenario: (OC/GC, baseline/colossus, crowding/no-crowding)
	'''
	# Getting paths for requisite data files
	# ocPath = os.path.join
	gcPath = ocPath.replace('Open', 'Globular')

	for ocFile, gcFile in zip(os.listdir(ocPath), os.listdir(gcPath)):
		
		# Ignoring wd file dsubdirs
		if ('.csv' in ocFile and '.csv' in gcFile) and (gcFile[0:4] == ocFile[0:4]): # and ('obs' in ocFile or 'rec' in ocFile):
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
			print('multi-colored corner plot made!')

# Looping through Open Cluster file path, no need to loop through OC and GC files 
# fucntion will pull data from and make plots for corresponding GC observing scenario files
for root, dirs, files in os.walk(os.path.join('clusters', 'OpenClusters'), topdown = True):
	for d in dirs:
		if 'data' in d:
			path_to_OCdata = os.path.join(root, d)
			makeMultiCornerPlot(path_to_OCdata)


# TODO: maybe include m10-m67 crossover plots (can't hurt)
# - Upload and run on Quest to save memory
# Add corner article citation to thesis refs





