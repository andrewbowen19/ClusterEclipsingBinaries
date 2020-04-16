'''
Testing script to create different kinds of plots
Ideally we will call this as an import module in other scripts

'''

import matplotlib.pyplot as plot
import corner
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.colors as colors

class PlotMe(object):
        '''
        Class that creates different plot types
        Each submethod is a different plot type
        '''

        def __init__(self, data):
                self.data = data

        def corner_plot(DataFrame, popType, name, contours = False):
                '''
                Function to make a corner plot with our different parameters -- using seaborns
                df - should be dataframe with values (all binary params) (each columnn is one element)
                popType - binary subpopulation type (All, Obs, Rec)
                name - name of scenario (Open/Globular/Long; Baseline/Colossus; Crowding/NoCrowding)
                contours - set as True if contours overlaid on plots desired, default: False
                '''
                # Plot with contours
                if contours == True:
                        print('Making corner plots with contours...')
                        df = DataFrame
                        f = corner.corner(df, labels = df.columns, bins = 20, plot_contours = True)
                        f.suptitle(popType + '-' + name + ' White Dwarfs', fontsize=24)
                        f.show()
                        f.savefig(f'./plots/whiteDwarfPairs/corner_plots/contours/{popType}-{name}-cornerPlot_contour_WDBinary.pdf')

                        print('Corner contour plots made!')

                # No contours
                elif contours == False:
                        print('Making corner plots...')
                        df = DataFrame
                        f = corner.corner(df, labels = df.columns, bins = 20,plot_contours = False)
                        f.suptitle(popType + '-' + name + ' White Dwarfs', fontsize = 24)
                        f.show()
                        f.savefig(f'./plots/whiteDwarfPairs/corner_plots/{popType}-{name}-cornerPlot_WDBinary.pdf')

                        print('Corner plots made!')

                print('On to the next!')
        



