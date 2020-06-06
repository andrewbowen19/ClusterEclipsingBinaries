# Script to read in 1st rows of data files in aaron's testing repo (OpSim RA, ID, Dec) - can use file to run on my laptop
# Goes in and gets 1st row of files (names) and then second row (data) from each data .csv file

# Trying to get all OpSim fields with this script - written to a file
import pandas as pd
import glob
import numpy as np
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from OpSim import OpSim
import csv
#######################################################################################################################

np.set_printoptions(threshold=np.nan)#Setting up so it prints all values to file

ops = {'ID': None,'RA': None,'Dec': None, 'Nobs' : None}

# Getting OpSim fields and coords from minion data file (OpSim Class from OpSim.py)
OpS = OpSim()
OpS.dbFile = '/Users/andrewbowen/CEB_Project/data/data_files/minion_1016_sqlite.db' #for the OpSim database   
OpS.getAllOpSimFields()

# Turning RA/Dec into coords (in degrees)
coords = SkyCoord(OpS.RA, OpS.Dec, unit=(u.degree, u.degree),frame='icrs')

# Wrapping coords
RAwrap = coords.ra.wrap_at(180.*u.degree).degree
Decwrap = coords.dec.wrap_at(180.*u.degree).degree

# Appending class instances to ops dict keys
ops['ID'] = np.append(ops['ID'], OpS.fieldID)
ops['RA'] = np.append(ops['RA'], OpS.RA)
ops['Dec'] = np.append(ops['Dec'], OpS.Dec)
ops['Nobs'] = np.append(ops['Nobs'], OpS.Nobs)

# turning dict into a dataframe
ops_df = pd.DataFrame(ops, columns = ['ID','RA','Dec', 'Nobs'])

# Writing relevant bits (ID, RA, DEC) to a csv file
ops_df.to_csv('/Users/andrewbowen/CEB_Project/data/data_files/all_opsim_data.csv', sep = ' ', index = False)

