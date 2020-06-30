# Script to create plots from output eblsst files
# These files are grabbed from different observing scenarios
# eblsst files from the analyseClusterLISA script in the flex directory
# Let's plot!
# info on VRO filters here: 
# https://community.lsst.org/t/lsst-filter-profiles/1463

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re


class eblsstPlotter(object):
    '''
    Class to create histograms for different binary parameters
    This histogram data is pulled from EBLSST output files
    Parameters:
    DataFrame: pandas df style object pulled from csv files
    clusterType: either 'GlobularClusters', 'OpenClusters', 'm10', 'm67'
    crowding: str, 'Crowd' or 'NoCrowd'
    strategy: str, 'base', 'crowd'
    The string inputs are case-sensitive
    '''
    def __init__(self, DataFrame, clusterType, strategy, crowding, param):
        self.DataFrame = DataFrame
        self.clusterType = clusterType
        self.strategy = strategy
        self.param = param  # binary parameter to be plotted
        # Filter plotting colors '<filter name>': '<matplotlib color>'
        self.filterColors = {'u': '#8B00FF', 'g': '#0000FF',
                             'r': '#00FF00', 'i': '#FFFF00',
                             'z': '#FF7F00', 'y': '#FF0000',
                             'ugrizy': '#000000'}
        # Dictionary including param keys (from file names)
        self.paramDict = {'lp': 'log-period', 'd': 'dist',
                          'm1': 'm1', 'q': 'q',
                          'mag': 'mag', 'e': 'e', 'r': 'r'}
        self.crowding = crowding

    def makeIndHist(self, data, filter, show_mean=False, save_fig=True):
        '''
        Generates histograms of binary paramsfor individual filters
        data: list-like containing binaries in each bin 
            (should be a 'hist' column from eblsst files)
        bins: bins for hist, use 'binsEdges' columnn from csv files
        '''
        print(f'Making Individual histogram for {self.clusterType}-{self.strategy}-{self.param}...')
        print(f'filter: {filter}')
        f, ax = plt.subplots(figsize=(8, 5))
        ax.hist(data, bins=self.DataFrame['binEdges'],
                histtype='step', color=self.filterColors[filter])
        ax.set_title(self.clusterType + '-' + self.strategy + self.crowding + '; ' + filter)
        ax.set_xlabel(self.param)

        if show_mean:
            ax.vline(np.mean(data), linestyle='--', label='Distribution mean')

        if save_fig:
            fileName = f'{self.clusterType}-{self.strategy}-{self.param}-{filter}-histInd.pdf'
            path_to_fig = os.path.join('.', 'plots', self.clusterType, self.strategy + self.crowding, 'ind', fileName)
            f.savefig(path_to_fig, format='pdf')

        print('Individual histogram complete!')

    def makeCombinedHist(self, filters, save_fig=True):
        '''
        Method to make histogram of binaries including filters
        data: just pass dataframe with filter columns to this
        filters: list of strings containinf filter keys
        '''
        multiData = self.DataFrame.drop(['binEdges', 'histAll', 'histObs', 'allhistRec'],
                                        axis=1)

        fig, ax = plt.subplots(figsize=(8, 5))
        # Looping thru filters
        for f in filters:
            ax.hist(multiData[f+'_histRec'], bins=self.DataFrame['binEdges'],
                    histtype='step', linestyle='--',
                    color=self.filterColors[f], label=f, alpha=0.7)
        ax.legend()
        ax.set_title(ax.set_title(self.clusterType + '-' + self.strategy + self.crowding + '; recovered'))
        ax.set_xlabel(self.paramDict[self.param])  # replacing parameter label from file name w param key (in dict above)
        ax.set_ylabel(r'$N_{bin}$')

        if save_fig:
            filterStr = ''.join(filters)
            fileName = f'{self.clusterType}-{self.strategy}-{self.param}-{filterStr}-histCombined.pdf'
            path_to_fig = os.path.join('.', 'plots', self.clusterType, self.strategy + self.crowding, 'combined', fileName)
            fig.savefig(path_to_fig)


