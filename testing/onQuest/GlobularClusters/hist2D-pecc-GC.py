# Andrew's plotting script for period-eccentricity tidal circularization plots
# Want scatter plot of eccentricity vs period colored by Nrec
# 2D histogram (CIERA REU Python sequence part 5) - Nrec and %rec in each ecc-p bin

# Globular Clusters Version -- ON Quest

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.cm as cm

# datAll = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/all-ecc-p.csv', header = 0, names = ['e', 'p'])
# datObs = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/obs-ecc-p.csv', header = 0, names = ['e', 'p'])
# datRec = pd.read_csv('/Users/andrewbowen/ceb_project/testing/GlobularClusters/pecc/rec-ecc-p.csv', header = 0, names = ['e', 'p'])

# Test 2D hist - 'all' dataframe
# plt.hist2d(datAll['p'], datAll['e'], bins = [10,10])
# plt.xlabel('period (days)')
# plt.ylabel('eccentricity')
# plt.show()

##################################################################################################################

def scatterPecc(df, plot_title):
    # Function to make a scatter plot of p-ecc as a check for our population stats
    p = df['p']
    ecc = df['e']
    
    f,ax = plt.subplots()
    ax.scatter(np.log10(p),ecc)
    ax.set_xlabel('period (log-days)')
    ax.set_ylabel('ecc')
    ax.set_title(plot_title)
    # Add in plot tile in for loop below - want these scatters for 3 subpops and all scenarios (O/G, B/C, C/N)
    f.savefig('./plots/' + plot_title + '-scatter-pecc.pdf')
    
############################################################################################


def PercentRecHist(recHist, obsHist, title):
	#Function to make percent Recovered histogram of period and eccentricity - going to put on quest scripts
	trueRecHist = np.divide(recHist,obsHist)
	print('Obs: ', obsHist, 'Rec: ', recHist)
	# Plotting percent rec hist
	f, ax = plt.subplots()
	
	im = ax.imshow(np.nan_to_num(trueRecHist), cmap = cm.Blues,aspect = (xmax - xmin)/(ymax - ymin),origin = 'lower', extent = [xmin,xmax,ymin, ymax])
	ax.set_xlabel('Period (log-days)')
	ax.set_ylabel('Eccentricity')
	ax.set_title(r'% recovered: ' + title)
	plt.show()
	plt.savefig('./plots/contour_plots/percentRecHist'+ title + '.pdf')
    




##################################################################################################################

# Function copied from Gaia_M67_SOLUTIONS (/CIERA_REU/PythonTutorialSequence/Part5/)
# Want to convert from plotting proper motion to ecc and p
# Can call this for each scenario (GC/OC, baseline/colossus, crowding/noCrowding)
def plotPecc(p, ecc,  xmin, xmax, Nx, ymin, ymax, Ny, plot_title, norm = None):
    f = plt.figure(figsize=(8, 8)) 
    gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1]) 
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])
    print('log-period: ', np.log10(p))
    #histograms
    pBins = np.linspace(xmin,xmax,Nx)
    eccBins = np.linspace(ymin,ymax,Ny)
    hx1D, x1D, im = ax1.hist(np.log10(p), bins=pBins, histtype='step', fill=False)
    hy1D, y1D, im = ax3.hist(ecc, bins=eccBins, histtype='step', fill=False, orientation="horizontal")

    #heatmap
    #a		x2.hexbin(r['p']*np.cos(r['dec']*np.pi/180.), r['ecc'], gridsize=30, cmap=cm.Blues, bins='log', extent=(xmin, xmax, ymin, ymax))
    h2D, x2D, y2D, im = ax2.hist2d(np.log10(p), ecc, bins=[Nx, Ny], \
                                    range=[[xmin, xmax], [ymin, ymax]], norm = norm, cmap = cm.Blues)
    ax1.set_xlim(xmin, xmax)
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    #ax2.set_xscale('log')
    ax3.set_ylim(ymin, ymax)
    ax2.set_xlabel(r'period (log-days)', fontsize=16)
    ax2.set_ylabel(r'binary eccentricity', fontsize=16)
    plt.setp(ax1.get_yticklabels()[0], visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels()[0], visible=False)
    plt.title(plot_title, fontsize = 22) # Can add titles to distinguish between plots
    f.subplots_adjust(hspace=0., wspace=0.)
    # plt.show()
    f.savefig('./plots/contour_plots/'+ plot_title + '-p-ecc-hist2D.pdf')
    return x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, (ax1, ax2, ax3)



# ########################################################################################################################
# ########################################################################################################################

# Setting up bins: will likely change
xmin, xmax, Nx = -1, 2, 100
ymin, ymax, Ny = 0, 1, 100

# Testing with local files - baseline GCs no crowding
# All
# x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, ax = plotPecc(datAll['p'], datAll['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'All',norm=mpl.colors.LogNorm())
# Obs
# x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, ax = plotPecc(datObs['p'], datObs['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'Obs',norm=mpl.colors.LogNorm())
# Rec
#x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, ax = plotPecc(datRec['p'], datRec['e'], xmin, xmax, Nx, ymin, ymax, Ny, 'Rec',norm=mpl.colors.LogNorm())
#print(np.max(h2D.flatten()))

# Could try looping thru all the paths to the all,obs,rec files for each scenario
import os
n = 0

# Walks through every file in directory, only pulls relevant csv files and makes plots for each scenario
for root, dirs, files in os.walk(".", topdown=True):
	for name in files:
		# Only reading in p-ecc output files
		if '-ecc-p.csv' in name:
			df = pd.read_csv(os.path.join(root,name), header = 0, names = ['e', 'p'])
			print(n)
			#print(df)
			n += 1
			x2D, y2D, h2D, x1D, hx1D, y1D, hy1D, ax = plotPecc(np.log10(df['p']),df['e'], xmin, xmax, Nx, ymin, ymax, Ny, name[0:7], norm=mpl.colors.LogNorm())
			
			scatterPecc(df, name[0:7])
			# Only want matching scenario obs-rec dataframes

			if ('obs' in name) and ('.csv' in name):
				obsHist2D = h2D
			elif 'rec' in name and '.csv' in name:
				recHist2D = h2D
			# Making 2D histgram of percent recoved
				PercentRecHist(recHist2D, obsHist2D, name[4:7])

		else:
			pass









