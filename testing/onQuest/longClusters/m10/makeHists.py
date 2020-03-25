# New plotting script to creat 2D histograms for different binary params (p, ecc, mass, inc, radius, lum, etc.)
# M10 Long clusters version -- will make plots for all scenarios and subpopulations 

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

import corner
from itertools import permutations

def make2dHist(x, y, xmin, xmax, Nx, ymin, ymax, Ny, popType, title, xlabel, ylabel, plotDir, xlog = False, norm = None):
	'''
	Function to make 2D histograms (Nrec/Obs/All), popType gives binary subpopulation (all, obs, rec)
	'''
	f = plt.figure(figsize=(8, 8)) 
	gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[2])
	ax3 = plt.subplot(gs[3])

	xbins = np.linspace(xmin,xmax,Nx)
	ybins = np.linspace(ymin,ymax,Ny)

	# Log option for plots with period
	if xlog == True:
		logx = np.log10(x)
		hx1D, x1D, im = ax1.hist(logx, bins = xbins, histtype = 'step', fill = False)
		hy1D, y1D, im = ax3.hist(y, bins = ybins, histtype = 'step', fill = False, orientation = "horizontal")

		h2D, x2D, y2D, im = ax2.hist2d(logx, y, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = colors.LogNorm(vmin = 0.1, vmax = 1), cmap = cm.Blues)

		# Setting xlabel if we're plotting in the log
		ax2.set_xlabel(f'log({xlabel})', fontsize=16)


	# Plotting without log
	else:
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


def PercentRecHist(xRec, yRec, xObs, yObs, xmin, xmax, Nx, ymin, ymax, Ny, popType, title, xlabel, ylabel, plotDir ,xlog = False, norm = None):	
	'''
	Fucntion that takes in raw xRec and 
	'''
	f = plt.figure(figsize=(8, 8)) 
	gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[2])
	ax3 = plt.subplot(gs[3])

	xbins = np.linspace(xmin,xmax,Nx)
	ybins = np.linspace(ymin,ymax,Ny)

	# xBins = np.linspace(xmin, xmax, Nx)
	# yBins = np.linspace(ymin, ymax, Ny)

	xRecHist, xRbins = np.histogram(xRec, bins = xbins)
	xObsHist, xObins = np.histogram(yRec, bins = xbins)

	yRecHist, yRbins = np.histogram(yRec, bins = ybins)
	yObsHist, yObins = np.histogram(yObs, bins = ybins)

	xPercentRec = np.nan_to_num(np.divide(xRecHist, xObsHist))
	yPercentRec = np.nan_to_num(np.divide(yRecHist, yObsHist))

	## print('Recovered Hists: ', xPercentRec, yPercentRec)

	hy1D, y1D, im = ax3.hist(yPercentRec, bins=yRbins, range = (xmin, xmax),histtype='step', fill=False, orientation="horizontal")	

		# Log option for plots with period
	if xlog == True:
		logxRec = np.log10(xRec)
		logxObs = np.log10(xObs)
		hx1D, x1D, im = ax1.hist(np.log10(xPercentRec), bins=xRbins, range = (xmin, xmax),histtype='step', fill=False)

		# Rec hist 2d
		h2D_rec, x2D_rec, y2D_rec, im = ax2.hist2d(logxRec, yRec, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)

		# Obs hist 2d
		h2D_obs, x2D_obs, y2D_obs, im = ax2.hist2d(logxObs, yObs, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)

		percentRecHist = np.divide(np.nan_to_num(h2D_rec), np.nan_to_num(h2D_obs))
		# Setting xlabel if we're plotting in the log
		ax2.set_xlabel(f'log({xlabel})', fontsize=16)

	# if we don't want the x-scale in the log (not plotting period)
	else:
		hx1D, x1D, im = ax1.hist(xPercentRec, bins=xRbins, range = (xmin, xmax),histtype='step', fill=False)
		# Rec hist 2d
		h2D_rec, x2D_rec, y2D_rec, im = ax2.hist2d(xRec, yRec, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)
		# Obs hist 2d
		h2D_obs, x2D_obs, y2D_obs, im = ax2.hist2d(xObs, yObs, bins=[Nx, Ny], \
			range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)

		percentRecHist = np.divide(h2D_rec,h2D_obs)
		ax2.set_xlabel(xlabel, fontsize=16)

		# Plotting actual 2d % rec hist - need to reanspose with .T method
	im = ax2.imshow(np.nan_to_num(percentRecHist), cmap = cm.Blues,
		aspect = (xmax - xmin)/(ymax - ymin), origin = 'lower', extent = [xmin,xmax,ymin, ymax], norm = colors.LogNorm(vmin = 0.1, vmax = 1))
	ax2.set_ylabel(ylabel, fontsize=16)
	plt.setp(ax1.get_yticklabels()[0], visible=False)
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.setp(ax3.get_yticklabels(), visible=False)
	plt.setp(ax3.get_xticklabels()[0], visible=False)
	plt.title(title, fontsize = 22) # Can add titles to distinguish between plots
	f.subplots_adjust(hspace=0., wspace=0.)
	# Saving fig
	f.savefig(f'./plots/hists2D/percentRec/{plotDir}/{popType}-{title}-{xlabel}-{ylabel}.pdf')

# Maybe put this in its own script for runtime purposes
def corner_plot(df, popType, name):
	'''
	Function to make a corner plot with our different parameters -- using seaborns
	df - should be dataframe with values (all binary params) (each columnn is one element)
	popType - binary subpopulation type (All, Obs, Rec)
	name - name of scenario (Open/Globular/Long; Baseline/Colossus; Crowding/NoCrowding)
	'''
	print('Making corner plots...')

	df['p'] = np.log10(df['p']) # Plotting p in the log
	f = corner.corner(df, labels = df.columns, bins = 50, show_contours = False)
	f.show()
	f.savefig(f'./plots/corner_plots/{popType}-{name}-cornerPlot.pdf')

	print('Corner plots made!')


### Corner test ###

# dat = pd.read_csv('./data/all-M10BN-histData.csv')
# dat = dat.drop('Unnamed: 0', axis =1) # idk where this random 'Unnamed: 0' column is coming from, probably when we write the files -- I think it's an index thing

# print('File read-in:  ',dat)

# corner_plot(dat, 'All', 'M10BN', None)



# ####################################################################################################################
# Keeping bin # at 100
xmin, xmax, Nx = 0, 2, 100
ymin, ymax, Ny = 0, 1, 100

# Setting up dictionary for labels and their corresponding axes limits; lists go by (xmin, xmax,  ymin, ymax)
limDict = {'p': [0.01,100], 'e':[0.,1.], 'm1':[0., 10.], 'm2':[0., 10.], 'i':[0.,100.], 'r1':[0.0, 10.], 'r2':[0.0, 10.]}

# Dictionary for which directory to save to
dirDict = {'CN': 'Col-NoCrowd', 'CC':'Col-Crowd', 'BC':'Base-Crowd', 'BN':'Base-NoCrowd'}


# looping thru directory files to find the correctly formatted csv files
for root, dirs, files in os.walk('.', topdown = True):

	haveObs = False
	haveRed = False
	for name in files:

		print('NAME: ',name)
		if '-histData.csv' in name:

			df = pd.read_csv(os.path.join(root,name), header = 0)
			df = df.drop('Unnamed: 0', axis =1)

			print('Data read in...')
			
			# Making corner plots for every scenario -- add to other makeHists scripts (M67, OCs and GCs)
			corner_plot(df, name[0:3], name[4:9])

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
				 name[0:3], name[4:9], pair[0], pair[1], dirDict[name[7:9]], xlog = xLog, norm = None)



			# Now making % recovered 2d histograms
			if 'obs' in name:
				df_obs = df
				haveObs = True
			
			elif 'rec' in name:
				df_rec = df
				haveRec = True
			# print('fooo foooo foooo: ',df_obs)

			if (haveObs and haveRec):	

				for pair in permutations(df_obs.columns,2): # using df here because they all will have the same columns
					# Do period plots in log
					if pair[0] == 'p':
						xLog = True
					else:
						xLog = False
					print(r'creating 2d % recHists...')

					PercentRecHist(df_rec[pair[0]], df_rec[pair[1]], df_obs[pair[0]], df_obs[pair[1]], limDict[pair[0]][0], limDict[pair[0]][1], Nx, limDict[pair[1]][0], limDict[pair[1]][1], Ny,
						name[0:3], name[4:9], pair[0], pair[1], dirDict[name[7:9]], xlog = xLog, norm = None)



print('All done!')


	

# PercentRecHist(xRec, yRec, xObs, yObs, xmin, xmax, Nx, ymin, ymax,
# Ny, popType, title, xlabel, ylabel, plotDir ,xlog = False, norm = None):	




# ADD log params 

# Corner plot - look into searborn (corner_plot)

# add magnitudes also











