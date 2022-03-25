# Corner plot script 
# Mean to plot for LISA WD binaries , magnitudes included in these plots
# This script should run for all 4 cluster types and observing scenarios. 
# Upload corner plots to box (both contours and no contour versions)
# Pulls from LISA files (may not be enough lines in obs and rec files, if so pass)

# Corner plotting documentation here: https://corner.readthedocs.io/en/latest/

import corner
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os

# Putting these functions this in its own script for runtime purposes
def corner_plot(DataFrame: pd.DataFrame, popType: str, name: str, contours: bool = False):
	'''
	Function to make a corner plot with our different parameters -- using seaborns
	parameters:
		df : pd.DataFrame; contain binary param values (all binary params) (each columnn is one element)
		popType : str; binary subpopulation type (All, Obs, Rec)
		name : str; name of scenario (Open/Globular/Long; Baseline/Colossus; Crowding/NoCrowding)
		contours : bool, default False; set as True if contours overlaid on plots desired, default: False
	'''
	# Plot with contours
	if contours == True:
		print('Making corner plots with contours...')
		df = DataFrame
		f = corner.corner(df, labels = df.columns, bins = 20, plot_contours = True)
		f.suptitle(popType + '-' + name + ' White Dwarfs', fontsize=24)
		f.show()
		f.savefig(f'./plots/lisaPlots/contours/{popType}-{name}-cornerPlot_contour_WDBinary.pdf')
		f.close()
		print('Corner contour plots made!')

	# No contours
	elif contours == False:
		print('Making corner plots...')
		df = DataFrame
		f = corner.corner(df, labels = df.columns, bins = 20, plot_contours = False)
		f.suptitle(popType + '-' + name + ' White Dwarfs', fontsize = 24)
		f.show()
		f.savefig(f'./plots/lisaPlots/{popType}-{name}-cornerPlot_WDBinary.pdf')
		f.close()
		print('Corner plots made!')

	print('On to the next!')



# ########################################################################################################

if __name__ == "__main__":
	# Looping through correct files in our trees and making plots
	for root, dirs, files in os.walk('./lisaFiles/', topdown = True):

		for name in files:

			print('ROOT, NAME: ', root, name)
			if 'WD-histDataLISA.csv' in name:

				dat = pd.read_csv(os.path.join(root,name), header = 0)
				# dat = dat.drop('Unnamed: 0', axis =1)
				# dat['p'] = np.log10(dat['p'])
				if len(dat) > 1:
					dat = dat.loc[np.where(dat['appMagMean_r'] != -999.0)] # Only want wds with good magnitudes
					dat.columns = ['p(days)', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$', 'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)','appMagMean_r']

					print('Data read in...')
					
					# Making corner plots for every scenario -- add to other makeHists scripts (M67, OCs and GCs)
					if len(dat) > len(dat.columns):
						# Using only relevant slices of names for OCs and GCs
						if ('GlobularClusters' in root) or ('OpenClusters' in root):
							corner_plot(dat, name[0:3], name[4:7], False) # No contours -- Change name slice to [4:7] for OCs and GCs (less characters in names)
							corner_plot(dat, name[0:3], name[4:7], True) # Making plots -- WITH contours

						# Not including long clusters for WD binaries -- thesis
						# If long clusters (name longer by 2 chars)
						elif ('m10' in root) or ('m67' in root):
							corner_plot(dat, name[0:3], name[4:9], False) # No contours -- Change name slice to [4:7] for OCs and GCs (less characters in names)
							corner_plot(dat, name[0:3], name[4:9], True) # Making plots -- WITH contours






