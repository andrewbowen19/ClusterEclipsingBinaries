#Script to read in 1st rows of data files in aaron's testing repo - can use file to run on my laptop
# Goes in and gets 1st row of files (names) and then second row (data) from each data .csv file


import pandas as pd
import glob
import numpy as np
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from OpSim import OpSim
import csv
# #####################################################################################################################
# Getting OpSim data
#np.set_printoptions(threshold=np.nan)

# IDdict['frac'] = np.append(IDdict['frac'],f)

ops = {'ID': None,'RA': None,'Dec': None}

OpS = OpSim()
OpS.dbFile = '/projects/p30137/abowen/CEB/minion_1016_sqlite.db' #for the OpSim database   
OpS.getAllOpSimFields()


coords = SkyCoord(OpS.RA, OpS.Dec, unit=(u.degree, u.degree),frame='icrs')

RAwrap = coords.ra.wrap_at(180.*u.degree).degree
Decwrap = coords.dec.wrap_at(180.*u.degree).degree

ops['ID'] = np.append(ops['ID'], np.array2string(OpS.fieldID))
ops['RA'] = np.append(ops['RA'], np.array2string(OpS.RA))
ops['Dec'] = np.append(ops['Dec'], np.array2string(OpS.Dec))

 

w = csv.writer(open("/projects/p30137/abowen/CEB/opsim_output.csv", "w"))
for key, val in ops.items():
	w.writerow([key, val])


# Importing nObs for LSST fields
new_path ='/projects/p30137/ageller/testing/EBLSST/fix_omega/output_files/'
allFiles = glob.glob(new_path + "/*.csv")

# Reading in 'example'
header_file = '/projects/p30137/ageller/testing/EBLSST/fix_omega/output_files/2140output_file.csv'
dat_obs = pd.read_csv(header_file, sep = ',', nrows = 0)#getting names from file

# reading in other frames
for file in allFiles:
	dat1 = pd.read_csv(file, sep = ',', header = 0, nrows = 1)#reads in: OpSimID,OpSimRA,OpSimDec,NstarsTRILEGAL,NOpSimObs_u,NOpSimObs_g,NOpSimObs_r,NOpSimObs_i,NOpSimObs_z,NOpSimObs_y
	dat_obs = dat_obs.append(dat1)


ID = dat_obs['OpSimID']
#dat_obs.set_index()
f = open('OpSim_data.txt', 'w')
f.write(dat_obs.to_string())
f.close()
