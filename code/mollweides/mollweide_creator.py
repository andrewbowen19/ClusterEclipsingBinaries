'''
Script to generate Number of observations, Nobs mollweides for colossus and baseline strategies

DO NOT LOSE THIS LINK:
Download site for different strategy db files
http://astro-lsst-01.astro.washington.edu:8080/?runId=16

Be sure to use colossus 2664 NOT 2665 which is baseline-style
'''


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
import sys

p = os.environ['PATH']
pv = os.path.join(os.getcwd(), 'vespa_update')
p2 = pv+':'+p
os.environ['PATH'] = p2

# Establishing path
# mypath = '/Users/andrewbowen/ClusterEclipsingBinaries/code/mollweides/'
mypath = os.getcwd()
sys.path.insert(0, mypath)

# OpSim viewing field diameter
fieldRad = Angle(1.75, unit=u.degree)


def plotMollweide(x, y, cVal, fieldSize, title):
    '''
    Function to generate mollweide projection plot
    cVal - values for colormap
    '''
    f, ax = plt.subplots(figsize=(8, 5),
                         subplot_kw={'projection': "mollweide"})
    ax.grid(True)
    ax.set_xlabel(r"$RA$", fontsize=16)
    ax.set_ylabel(r"$Dec$", fontsize=16)
    im = ax.scatter(x, y, s=fieldSize, c=cVal,
                    cmap='viridis', alpha=1, vmin=0, vmax=1000)

    ax.set_title(title, fontsize=24)
    cbar = f.colorbar(im)
    cbar.ax.set_ylabel(r'Nobs', fontsize = 16, rotation=270, labelpad=8)

    f.savefig(f'./plots/{title}-Nobs-mollweide.pdf', format='pdf', bbox_inches='tight')
    plt.show()


def wrapCoords(coords):
    '''
    Code from Aaron's script
    Function to reformat coordinates (wrapping)
    coords - astropy SkyCoord objects containing RA and Dec coords
    '''
    raGal = coords.icrs.ra.wrap_at(180.*u.degree).degree
    decGal = coords.icrs.dec.wrap_at(180.*u.degree).degree
    lGal = coords.galactic.l.wrap_at(180.*u.degree).degree
    bGal = coords.galactic.b.wrap_at(180.*u.degree).degree

    return raGal, decGal, lGal, bGal


def createMollweide(filename, strategy, create_csv=True):
    '''
    Function to create mollweide Nobs plots for different OpSim strategies
    '''
    print('Creating mollweides for: ', strategy)

    # Getting OpSim field data and Nobs for each field
    OpS = OpSim()
    OpS.dbFile = filename
    print(OpS.dbFile)
    OpS.getAllOpSimFields()

    # Setting up and wrapping OpSim coords
    coords = SkyCoord(OpS.RA, OpS.Dec, unit=(u.degree, u.degree), frame='icrs')
    raGal, decGal, lGal, bGal = wrapCoords(coords)

    xx = np.where(OpS.Nobs > 0)
    ra = raGal[xx]*np.pi/180.
    dec = decGal[xx]*np.pi/180.
    goodObs = OpS.Nobs[xx]

    # Making the mollweide
    plotMollweide(ra, dec, goodObs, fieldRad, strategy)

    # Writing Nobs and coordinate data to csv output files
    if create_csv:
        df = pd.DataFrame({'OpSim ID': OpS.fieldID,
                           'OpSim RA': OpS.RA, 'OpSim Dec': OpS.Dec,
                           'Nobs': OpS.Nobs})
        df.to_csv('./output/' + strategy + '-OpSim-FieldData.csv', index=False)
        print('CSV file created!')


# Creating mollweides: from OpSim db files - only need to run once!
# createMollweide('colossus_2664.db', 'colossus', create_csv=True)
# createMollweide('baseline2018a.db', 'baseline', create_csv=True)
# createMollweide('kraken_2026.db', 'kraken', create_csv=True)


# New function to create mollweides
def csvMollweide(filename, strategy):
    '''
    Function to create mollweides from csv files creaeted by above function
    Faster than instantiating OpSim for each file
    '''
    df = pd.read_csv(f'./output/{filename}')

    # Pulling only good Nobs fields and setting up coordinates for plotting
    dat = df.loc[df['Nobs'] > 0.0]
    coords = SkyCoord(dat['OpSim RA'], dat['OpSim Dec'],
                      unit=(u.degree, u.degree), frame='icrs')

    # Wrapping coords
    raGal, decGal, lGal, bGal = wrapCoords(coords)

    # Making the mollweides w/ wrapped coords
    plotMollweide(raGal*np.pi/180., decGal*np.pi/180.,
                  dat['Nobs'], fieldRad, strategy)


csvMollweide('colossus-OpSim-FieldData.csv', 'colossus')
csvMollweide('baseline-OpSim-FieldData.csv', 'baseline')
csvMollweide('kraken-OpSim-FieldData.csv', 'kraken')
