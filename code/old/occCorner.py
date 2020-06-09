'''
Quick & dirty corner script for colCrowd scenario
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
		f1 = corner.corner(DataFrame1, color = 'r', labels = DataFrame1.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, range = ranges, hist2d_kwargs={'quiet':True})#, scale_hist = True, hist2d_kwargs = {'plot_density':True})  # , scale_hist = True)
		# Plotting secodn dataframe - GCs
		corner.corner(DataFrame2, color = 'b', labels = DataFrame2.columns, label_kwargs={"fontsize":18},
							bins = 20, plot_contours = True, title_kwargs={"fontsize": 28}, fig = f1, range = ranges, hist2d_kwargs={'quiet':True})# hist2d_kwargs = {'plot_density':True})
		f1.suptitle(popType + '-' + name, fontsize=24)
		# plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/', f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner contour plots made!')

gcPath = os.path.join('clusters', 'GlobularClusters', 'colCrowd', 'data', 'rec-GCC-histData.csv')
ocPath = os.path.join('clusters', 'OpenClusters', 'colCrowd', 'data', 'rec-OCC-histData.csv')

gcDat = pd.read_csv(gcPath, header=0)
ocDat = pd.read_csv(ocPath, header=0)

corner_plot(ocDat, gcDat, 'rec', 'OCC-GCC', True)







