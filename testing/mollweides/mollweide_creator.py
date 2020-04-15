# Script to generate Nobs mollweides for colossus and baseline strategies
# Real version, one in 'code' directory is an imposter

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

def createMollweide(filename, strategy, create_csv=True):
    '''Function to create mollweide Nobs plots for different OpSim strategies'''
    print('Creating mollweides for ', strategy) 

    # Getting OpSim field data and Nobs for each field
    OpS = OpSim()
    OpS.dbFile = filename
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

    f, ax = plt.subplots(figsize = (8,5), subplot_kw={'projection': "mollweide"})
    ax.grid(True)
    ax.set_xlabel(r"$RA$", fontsize=16)
    ax.set_ylabel(r"$Dec$", fontsize=16)
    xx = np.where(OpS.Nobs > 0)
    ax.scatter(raGal[xx]*np.pi/180., decGal[xx]*np.pi/180., s=fieldRad, c=OpS.Nobs[xx], 
               cmap='tab10', alpha = 1, vmin=0, vmax=1000)
    ax.set_title(strategy, fontsize=16)
    f.savefig(f'./plots/{strategy}-Nobs-mollweide.pdf')
    plt.show()

    if create_csv:
        dat = pd.DataFrame({'OpSim ID' : OpS.fieldID,
                            'OpSim RA' : OpS.RA, 'OpSim Dec' : OpS.Dec, 'Nobs' : OpS.Nobs})
        dat.to_csv('./output/' + strategy + '-OpSim-FieldData.csv', index=False)
        print('CSV file created!')


# Creating mollweides
createMollweide('baseline2018a.db', 'baseline', create_csv=True)
createMollweide('colossus_2665.db', 'colossus', create_csv=True)

 