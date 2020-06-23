'''
Corner plotting script for plots with multiple colors for OCs/GCs
Loops through file path to create corner plots for each observing scenario
Updated for flex file paths (even for all observing scenarios)
Will create corner plots for all binary sample (i.e. not just WDs)
Corner plotting documentation: https://corner.readthedocs.io/en/latest/
Written a a function, no need for object in file tree loop
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


# ##########################################################################################################

def corner_plot(DataFrame1, DataFrame2, popType, name, contours=False):
	'''
	Function to make a corner plot with our different binary parameters
	Uses corner object to generate plots
	Creates corner plot for two different cluster sub-pops (OC vs GC)
	DataFrame1 - should be dataframe with values for all binary params
	DataFrame2 - same structure as df1, but for GCs
	popType - binary subpopulation type (All, Obs, Rec)
	name - name of observing scenario (OC/GC; Base/Col; Crowding/NoCrowd)
	contours (bool) - create plot w/ or w.out contours
	'''

	# Plot with Contours
	if contours is True:
		print('Making corner plots with contours...')
		# df1 = DataFrame1
		# df2 = DataFrame2

		# Setting up better range values for either OCs or GCs
		ranges = [(-0.5, 2.5), (0., 15.), (0, 15.),
				  (0., 8.), (0., 8.), (0., 1.),
				  (45., 135.), (15., 24.)]

		# Plotting first dataframe - OCs
		f1 = corner.corner(DataFrame1, color='r',
						   labels=DataFrame1.columns,
						   label_kwargs={"fontsize": 18},
						   bins=20, plot_contours=True,
						   title_kwargs={"fontsize": 28}, range=ranges,
						   hist2d_kwargs={'quiet': True})
		# Plotting secodn dataframe - GCs
		corner.corner(DataFrame2, color='b',
					  labels=DataFrame2.columns,
					  label_kwargs={"fontsize": 18}, bins=20,
					  plot_contours=True, title_kwargs={"fontsize": 28},
					  fig=f1, range=ranges, hist2d_kwargs={'quiet': True})
		f1.suptitle(popType + '-' + name, fontsize=24)
		# plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/',
								f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner contour plots made!')

	# No contours
	elif contours is False:
		print('Making corner plots...')

		# Plot 1
		f1 = corner.corner(DataFrame1, color='r', labels=DataFrame1.columns,
						   label_kwargs={"fontsize": 18}, bins=20,
						   plot_contours=False,
						   no_fill_contours=True,
						   plot_density=True, title_kwargs={"fontsize": 28},
						   hist_kwargs={"density": True})

		# Plot 2
		corner.corner(DataFrame2, color='b', labels=DataFrame2.columns,
					  label_kwargs={"fontsize": 18}, bins=20,
					  plot_contours=False,
					  no_fill_contours=True,
					  plot_density=True, title_kwargs={"fontsize": 28},
					  fig=f1, hist_kwargs={"density": True})
		f1.suptitle(popType + '-' + name, fontsize=24)
		# plt.show()
		f1.savefig(os.path.join('./plots/corner_plots/multi/',
				   f'{popType}-{name}-cornerPlot-contour.pdf'))
		plt.close()
		print('Corner plots made!')

	print('On to the next!')


# Looping through Open Cluster file path
# Function make plots for corresponding OC/GC files
for root, dirs, files in os.walk(os.path.join('.', 'clusters', 'OpenClusters'), topdown=True):
	for d in dirs:
		if 'data' in d:
			# Generating paths to files
			path1 = os.path.join(root, d)
			path2 = path1.replace('Open', 'Globular')

			for file1, file2 in zip(os.listdir(path1), os.listdir(path2)):
				print('FILES', file1, file2)

				# Ignoring wd file dsubdirs
				if ('.csv' in file1 and '.csv' in file2) and (file1[0:4] == file2[0:4]):

					# Reading in binary data for all 3 sub-pops
					dat1 = pd.read_csv(os.path.join(path1, file1), header=0)
					dat2 = pd.read_csv(os.path.join(path2, file2), header=0)

					# Converting period to log-days
					dat1['p'] = np.log10(dat1['p'])
					dat2['p'] = np.log10(dat2['p'])
					dat1.columns = ['log-p', 'm1 $(M_{\odot})$',
									'm2 $(M_{\odot})$',
									'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$',
									'e', 'i (deg)', 'App Mag Mean r']
					dat2.columns = ['log-p', 'm1 $(M_{\odot})$',
									'm2 $(M_{\odot})$',
									'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$',
									'e', 'i (deg)', 'App Mag Mean r']

					# Making multi corner plots - OCs in red, GCs in blue
					corner_plot(dat1, dat2, file1[0:3],
								file1[4:7] + '-' + file2[4:7], False)


# TODO: maybe include m10-m67 crossover plots (can't hurt)
# - Upload and run on Quest to save memory
# Add corner article citation to thesis refs
