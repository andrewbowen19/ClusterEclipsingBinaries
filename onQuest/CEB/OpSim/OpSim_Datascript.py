# Script to read in 1st rows of data files in aaron's testing repo - can use file to run on my laptop


import pandas as pd
import glob


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


f = open('OpSim_data.txt', 'w')
f.write(dat_obs.to_string())
f.close()
