# New plotting script to creat 2D histograms for different binary params (p, ecc, mass, inc, radius, lum, etc.)
# Quest Open Clusters Version


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
import os

#dat = pd.read_csv('./data/all-M10BN-histData.csv')
#print(dat)
print('########################################')


def make2dHist(x, y,xmin, xmax, Nx, ymin, ymax, Ny, popType, title, xlabel, ylabel, plotDir ,xlog = False, norm = None):
	'''Function to make 2D histograms (Nrec/Obs/All), popType gives binary subpopulation (all, obs, rec)'''
	f = plt.figure(figsize=(8, 8)) 
	gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[2])
	ax3 = plt.subplot(gs[3])

	xbins = np.linspace(xmin,xmax,Nx)
	ybins = np.linspace(ymin,ymax,Ny)
	# print('xbins: ', xbins)

	# Log option for plots with period
	if xlog == True:
		logx = np.log10(x)
		hx1D, x1D, im = ax1.hist(logx, bins = xbins, histtype = 'step', fill = False)
		hy1D, y1D, im = ax3.hist(y, bins = ybins, histtype = 'step', fill = False, orientation = "horizontal")

		h2D, x2D, y2D, im = ax2.hist2d(logx, y, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)

		# Setting xlabel if we're plotting in the log
		ax2.set_xlabel(f'log({xlabel})', fontsize=16)


	# Plotting without log
	elif xlog == False:
		hx1D, x1D, im = ax1.hist(x, bins = xbins, histtype = 'step', fill = False)
		hy1D, y1D, im = ax3.hist(y, bins = ybins, histtype = 'step', fill = False, orientation = "horizontal")

		h2D, x2D, y2D, im = ax2.hist2d(x, y, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)
		ax2.set_xlabel(xlabel, fontsize=16)

	# Adjusting and setting labels
	ax1.set_xlim(xmin, xmax)
	ax2.set_xlim(xmin, xmax)
	ax2.set_ylim(ymin, ymax)
	ax3.set_ylim(ymin, ymax)
	
	ax2.set_ylabel(ylabel, fontsize=16)
	plt.setp(ax1.get_yticklabels()[0], visible=False)
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.setp(ax3.get_yticklabels(), visible=False)
	plt.setp(ax3.get_xticklabels()[0], visible=False)
	plt.title(title, fontsize = 22) # Can add titles to distinguish between plots
	f.subplots_adjust(hspace=0., wspace=0.)
		# saving figure
	f.savefig(f'./plots/hists2D/{plotDir}/{popType}-{title}-{xlabel}-{ylabel}.pdf')


	return x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, (ax1, ax2, ax3)


def PercentRecHist():	
	print('foo')


# ####################################################################################################################
# Keeping bin # at 100
xmin, xmax, Nx = 0, 2, 100
ymin, ymax, Ny = 0, 1, 100

# Setting up dictionary for labels and their corresponding axes limits; lists go by (xmin, xmax,  ymin, ymax)
limDict = {'p': [0.01,100], 'e':[0.,1.], 'm1':[0., 10.], 'm2':[0., 10.], 'i':[0.,100.], 'r1':[0.0, 10.], 'r2':[0.0, 10.]}

# Dictionary for which directory to save to
dirDict = {'CN': 'Col-NoCrowd', 'CC':'Col-Crowd', 'BC':'Base-Crowd', 'BN':'Base-NoCrowd'}
# make2dHist(dat['p'], dat['m1'], xmin, xmax, Nx, ymin, ymax, Ny, 'all', 'M10BN', 'period', 'mass1', True, None)

from itertools import permutations

# looping thru directory files to find the correctly formatted csv files
for root, dirs, files in os.walk('.', topdown = True):


	for name in files:
		print('NAME: ',name)
		if '-histData.csv' in name:

			df = pd.read_csv(os.path.join(root,name),header = 0)
			df = df.drop('Unnamed: 0', axis =1)

			print('Data read in...')
			print('making plots...')

			# Making plots for each pair of params
			for pair in permutations(df.columns,2):
				# make2dHist(x, y,xmin, xmax, Nx, ymin, ymax, Ny, popType, title, xlabel, ylabel, xlog = False, norm = None)
				
				# Do period plots in log
				if pair[0] == 'p':
					xLog = True
				else:
					xLog = False
				print('creating 2d Hists...')

				# making raw # plots -- NOT % rec
				make2dHist(df[pair[0]], df[pair[1]], limDict[pair[0]][0], limDict[pair[0]][1], Nx, limDict[pair[1]][0], limDict[pair[1]][1], Ny,
				 name[0:3], name[4:7], pair[0], pair[1], dirDict[name[5:7]], xlog = xLog, norm = None)




print('All done!')
	
	

















