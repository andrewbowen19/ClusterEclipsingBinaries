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

class multiCorner(object):
	'''
	class that propduces and saves multi-colored corner plots
	Different colors correspond to different cluster types: OC/GC, M10/M67, etc.
	path - should be path to either OC or GC files, object pulls corresponding files from other cluster type
	'''

	def __init__(self, path):
		self.path = path
		# self.DataFrame1 = DataFrame1
		# self.DataFrame2 = DataFrame2
		# self.popType = popType
		# self.name = name
		# self.path2 = None

	def corner_plot(self, DataFrame1, DataFrame2, popType, name, contours = False):
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

		# No contours
		elif contours == False:
			print('Making corner plots...')
			df1 = DataFrame1
			df2 = DataFrame2
			f1 = corner.corner(DataFrame1, color = 'r', labels = DataFrame1.columns, label_kwargs={"fontsize":18},
								bins = 20, plot_contours = False, title_kwargs={"fontsize": 28})# , scale_hist = True)
			corner.corner(DataFrame2, color = 'b', labels = DataFrame2.columns, label_kwargs={"fontsize":18},
								bins = 20, plot_contours = False, title_kwargs={"fontsize": 28}, scale_hist = True, fig = f1)
			f1.suptitle(popType + '-' + name, fontsize=24)
			# plt.show()
			f1.savefig(os.path.join('./plots/corner_plots/multi/', f'{popType}-{name}-cornerPlot-contour.pdf'))
			plt.close()
			print('Corner plots made!')

		print('On to the next!')


	def makeMultiCornerPlot(self, path):
		'''
		Function to generate multi-colored corner plot
		Creates plot for an individual observing scenario: (OC/GC, baseline/colossus, crowding/no-crowding)
		'''
		# Getting paths for requisite data files: 
		# input path/cluster type corresponds to opposite cluster type for path2
		self.path2 = self.path
		if 'OpenClusters' in self.path:
			self.path2 = self.path.replace('Open', 'Globular')

		elif 'GlobularClusters' in self.path:
			self.path2 = self.path.replace('Globuular', 'Open')

		elif 'm10' in self.path:
			self.path2 = self.path.replace('m10', 'm67')

		elif 'm67' in self.path:
			self.path2 = self.path.replace('m67', 'm10')


		for file1, file2 in zip(os.listdir(self.path), os.listdir(self.path2)):
			
			# Ignoring wd file dsubdirs
			if ('.csv' in file1 and '.csv' in file2) and (file1[0:4] == file2[0:4]):  # and ('obs' in ocFile or 'rec' in ocFile):
				print('File name:', file1, file2)

				# Reading in binary data for all 3 sub-pops
				self.ocDat = pd.read_csv(os.path.join(self.path, file1), header = 0)
				self.gcDat = pd.read_csv(os.path.join(self.path2, file2), header = 0)
				print(len(self.ocDat), len(self.gcDat)) # checking to make sure file lenghts are in the same ballpark

				# Changing columns/period to log
				self.ocDat['p'] = np.log10(self.ocDat['p']) # converting period values from days to log-days
				self.gcDat['p'] = np.log10(self.gcDat['p']) # converting period values from days to log-days
				self.ocDat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$',
												'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r']
				self.gcDat.columns = ['log-p', 'm1 $(M_{\odot})$', 'm2 $(M_{\odot})$',
												'r1 $(R_{\odot})$', 'r2 $(R_{\odot})$', 'e', 'i (deg)', 'App Mag Mean r'] 
				# print('Binary Pop data for OCs, GCs: ', ocDat, '\n', gcDat)

				# Making multi corner plots - OCs in red, GCs in blue
				self.corner_plot(self.ocDat, self.gcDat, file1[0:3], file1[4:7] + '-' + file2[4:7], True)
				print('multi-colored corner plot made!')


# Test to see why OCC situation is being off
# new_path = './clusters/GlobularClusters/colCrowd/data/'
# mc_gcc = multiCorner(new_path)
# mc_gcc.makeMultiCornerPlot(new_path)

# Looping through Open Cluster file path, no need to loop through OC and GC files 
# fucntion will pull data from and make plots for corresponding GC observing scenario files
for root, dirs, files in os.walk(os.path.join('clusters', 'OpenClusters'), topdown = True):
	for d in dirs:
		print(d)
		if 'data' in d:
			print(root)
			path_to_OCdata = os.path.join(root, d)
			
			mc = multiCorner(path_to_OCdata)
			mc.makeMultiCornerPlot(path_to_OCdata)


# TODO: maybe include m10-m67 crossover plots (can't hurt)
# - Upload and run on Quest to save memory
# Add corner article citation to thesis refs





