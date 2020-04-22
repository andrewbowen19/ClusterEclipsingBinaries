# Script to generate Nobs mollweides for colossus and baseline strategies
# Real version, one in 'code' directory is an imposter

# DO NOT LOSE THIS LINK:
# Download site for different strategy db files
# http://astro-lsst-01.astro.washington.edu:8080/?runId=16

# Be sure to use colossus 2664 NOT 2665 which is baseline-style

import pandas as pd
import numpy as np
import os
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from OpSim import OpSim
from TRILEGAL import TRILEGAL
from matplotlib import pyplot as plt
from vespa_update import trilegal
import os
p = os.environ['PATH']
pv = os.path.join(os.getcwd(),'vespa_update')
p2 = pv+':'+p
os.environ['PATH'] = p2

import sys

mypath = '/Users/andrewbowen/ClusterEclipsingBinaries/testing/mollweides/'
sys.path.insert(0, mypath)


def plotMollweide(x, y, cVal, fieldSize, title):
    '''
    Function to generate mollweide projection plot
    cVal - values for colormap
    '''
    f, ax = plt.subplots(figsize = (8,5), subplot_kw={'projection': "mollweide"})
    ax.grid(True)
    ax.set_xlabel(r"$RA$", fontsize=16)
    ax.set_ylabel(r"$Dec$", fontsize=16)
    im = ax.scatter(x, y, s=fieldSize, c=cVal, 
               cmap='tab10', alpha = 1, vmin=0, vmax=1000)
    ax.set_title(title, fontsize=16)
    cbar = f.colorbar(im)
    cbar.ax.set_ylabel(r'Nobs', rotation=270, labelpad=8)

    f.savefig(f'./plots/{title}-Nobs-mollweide.pdf')
    plt.show()

def createMollweide(filename, strategy, create_csv=True):
    '''Function to create mollweide Nobs plots for different OpSim strategies'''
    print('Creating mollweides for: ', strategy) 

    # Getting OpSim field data and Nobs for each field
    OpS = OpSim()
    OpS.dbFile = filename
    print(OpS.dbFile)
    OpS.getAllOpSimFields()

    # OpSim Coords
    coords = SkyCoord(OpS.RA, OpS.Dec, unit=(u.degree, u.degree), frame='icrs')
    # OpSim field radius
    fieldRad = Angle(1.75, unit=u.degree)

    # Code from Aaron's script
    raGal = coords.icrs.ra.wrap_at(180.*u.degree).degree
    decGal = coords.icrs.dec.wrap_at(180.*u.degree).degree
    lGal = coords.galactic.l.wrap_at(180.*u.degree).degree
    bGal = coords.galactic.b.wrap_at(180.*u.degree).degree

    xx = np.where(OpS.Nobs > 0)
    ra = raGal[xx]*np.pi/180.
    dec = decGal[xx]*np.pi/180.
    goodObs = OpS.Nobs[xx]

    plotMollweide(ra, dec, goodObs, fieldRad, strategy)
   

    if create_csv:
        dat = pd.DataFrame({'OpSim ID' : OpS.fieldID,
                            'OpSim RA' : OpS.RA, 'OpSim Dec' : OpS.Dec, 'Nobs' : OpS.Nobs})
        dat.to_csv('./output/' + strategy + '-OpSim-FieldData.csv', index=False)
        print('CSV file created!')


# Creating mollweides
# createMollweide('colossus_2664.db', 'colossus', create_csv=True) # colossus first bc baseline works
# createMollweide('baseline2018a.db', 'baseline', create_csv=True)


# New function to create mollweides
def csvMollweide(filename, strategy):
    '''
    Function to create mollweides from csv files creaeted by above function
    '''
    df = pd.read_csv(f'./output/{filename}')

    # Pulling only good Nobs fields
    dat = df.loc[df['Nobs'] > 0.0]
    print(len(dat))

    coords = SkyCoord(dat['OpSim RA'], dat['OpSim Dec'], unit=(u.degree, u.degree), frame='icrs')
    # OpSim field radius
    fieldRad = Angle(1.75, unit=u.degree)

    # Code from Aaron's script - USED ABOVE (DRY)
    raGal = coords.icrs.ra.wrap_at(180.*u.degree).degree
    decGal = coords.icrs.dec.wrap_at(180.*u.degree).degree
    lGal = coords.galactic.l.wrap_at(180.*u.degree).degree
    bGal = coords.galactic.b.wrap_at(180.*u.degree).degree

    # Making the mollweides
    plotMollweide(raGal*np.pi/180., decGal*np.pi/180., dat['Nobs'], fieldRad, strategy)

csvMollweide('colossus-OpSim-FieldData.csv', 'colossus')
csvMollweide('baseline-OpSim-FieldData.csv', 'baseline')






