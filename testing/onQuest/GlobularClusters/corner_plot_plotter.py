# Corner plot script 
# Globular Cluster version, need to add to other file trees (M67, OC, GC)

import corner
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
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
		f = corner.corner(df, labels = df.columns, bins = 20, plot_contours = True)
		f.suptitle(popType + '-' + name, fontsize=24)
		f.show()
		f.savefig(f'./plots/corner_plots/contours/{popType}-{name}-cornerPlot-contour.pdf')

		print('Corner contour plots made!')

	# No contours
	elif contours == False:
		print('Making corner plots...')
		df = DataFrame
		f = corner.corner(df, labels = df.columns, bins = 20,plot_contours = False)
		f.suptitle(popType + '-' + name, fontsize = 24)
		f.show()
		f.savefig(f'./plots/corner_plots/{popType}-{name}-cornerPlot.pdf')

		print('Corner plots made!')

	print('On to the next!')

	

# ########################################################################################################

# Looping through correct files in our trees and making plots
for root, dirs, files in os.walk('.', topdown = True):

	for name in files:

		print('NAME: ',name)
		if '-histData.csv' in name:

			# Reading in an renaming dataframe columns for plotting - doing once and outside of our function call
			dat = pd.read_csv(os.path.join(root,name), header = 0)
			dat = dat.drop('Unnamed: 0', axis =1)
			dat['p'] = np.log10(dat['p'])
			dat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$', 'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)']


			print('Data read in...')
			
			# Making corner plots for every scenario -- add to other makeHists scripts (M67, OCs and GCs)
			corner_plot(dat, name[0:3], name[4:7], True) # WITH contours -- Change name slice to [4:7] for OCs and GCs (less characters)
			corner_plot(dat, name[0:3], name[4:7], False) # NO contours
			

print('All done!')





