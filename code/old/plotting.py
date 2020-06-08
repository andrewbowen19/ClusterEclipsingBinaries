"""

 New plotting script bc ithe one from last summer is bad and makes me angry

 """

######################################## script from final_plots.py - BAD BAD BAD ###################################
#																													#
#																													#
#																													#	
# # Importing needed libraries																						#
# import pandas as pd 																								#
# import numpy as np 																								#
# import matplotlib 																								#
# matplotlib.use('agg') 																							#
# import matplotlib.pyplot as plt 																					#
# # mattplotlib.use('agg')																							#
# from astropy import units as u 																					#
# from astropy.coordinates import SkyCoord, Distance 																#
#																													#
# # Reading in all data files at once 																				#
# import glob																										#
# # path = '/Users/andrewbowen/BO_EB-LSST/code/Opsim'																#
# path ='/projects/p30137/ageller/testing/EBLSST/fix_omega/output_files' # use your path							#
# allFiles = glob.glob(path + "/*.csv")																				#
# frame = pd.DataFrame()																							#
# frame1 = pd.DataFrame()																							#
# frame2 = pd.DataFrame()																							#
# list_1 = []																										#
# list_2 = []																										#
# IDdict = {'RA':np.array([]), 'Dec':np.array([]), 'frac':np.array([])} #For mollewiede plot later					#
#																													# 
# for file_ in allFiles:																							#
#     #print(file_)																									#
#     dat1 = pd.read_csv(file_, sep = ',', header=0, nrows = 1)														#
#     dat2 = pd.read_csv(file_, sep = ',', header=2)																#
#     dat2['RA'] = dat1['OpSimRA'][0]																				#
#     dat2['Dec'] = dat1['OpSimDec'][0]																				#
#     list_1.append(dat1)																							#
#     #list_2.append(dat2)																							#
#    																												#
#     if 'p' in dat2.keys():																						#
#         #print(dat2['p'].values[0])																				#
#         if dat2['p'].values[0] != -1:																				#	
#             list_2.append(dat2)																					#
#             IDdict['RA'] = np.append(IDdict['RA'],dat1['OpSimRA'])												#
#             IDdict['Dec'] = np.append(IDdict['Dec'], dat1['OpSimDec'])											#
#             PeriodIn = dat2['p']																					#
#             PeriodOut = dat2['LSM_PERIOD']																		#
#             pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values										#
#             f = 0																									#
#																													# 
#             recovered_p = dat2['LSM_PERIOD'].loc[pdiff<0.1]														#
#             all_p = dat2['LSM_PERIOD']																			#
# #     Putting recovered period into 'frac' lists																	#
#             f = recovered_p.shape[0]/all_p.shape[0]																#
# 																													#
# 																													#	
#             IDdict['frac'] = np.append(IDdict['frac'],f)															#
# frame1 = pd.concat(list_1) #RA, Dec and id stuff																	#
# frame2 = pd.concat(list_2) #Actual DataFrame 																		#
# ###################################################### what was I thinking???? ####################################

# Why didn't you comment more you dummy?


import pandas as pd 																								
import numpy as np 																								
import matplotlib 																								
# matplotlib.use('agg') 																						
import matplotlib.pyplot as plt 																					
# mattplotlib.use('agg')																							
from astropy import units as u 																					
from astropy.coordinates import SkyCoord, Distance 	

import glob																										
# path = '/Users/andrewbowen/BO_EB-LSST/code/Opsim'																
path ='/Users/andrewbowen/ceb_project/data/test_files/' # use your path							
allFiles = glob.glob(path + "/*.csv")																				
# frame = pd.DataFrame()																							
# frame1 = pd.DataFrame()																							
# frame2 = pd.DataFrame()																							
# list_1 = []																										
# list_2 = []																										
# IDdict = {'RA':np.array([]), 'Dec':np.array([]), 'frac':np.array([])} #For mollewiede plot later	

# for file_ in allFiles:																							
#     #print(file_)																									
#     dat1 = pd.read_csv(file_, sep = ',', header=0, nrows = 1)#reading in files - only header I guess											
#     dat2 = pd.read_csv(file_, sep = ',', header=2)#reading in 2nd header
#     dat2['RA'] = dat1['OpSimRA'][0]																				
#     dat2['Dec'] = dat1['OpSimDec'][0]																				
#     list_1.append(dat1)																							
#     list_2.append(dat2)																				

#     if 'p' in dat2.keys():																						
#         #print(dat2['p'].values[0])																				
#         if dat2['p'].values[0] != -1:																				
#             list_2.append(dat2)																					

# Names from OpSim output files
opsim_names = ['OpSimID','OpSimRA','OpSimDec','NstarsTRILEGAL','NOpSimObs_u','NOpSimObs_g','NOpSimObs_r',\
'NOpSimObs_i','NOpSimObs_z','NOpSimObs_y']

bin_names = ['p','m1','m2','r1','r2','e','i','d','nobs','Av',\
'[M/H]','appMagMean','maxDeltaMag','mag_failure','incl_failure','period_failure','radius_failure',\
'u_LSS_PERIOD','g_LSS_PERIOD','r_LSS_PERIOD','i_LSS_PERIOD','z_LSS_PERIOD','y_LSS_PERIOD','LSM_PERIOD']


# Read in 1 file, grab headers, then read in every file and concat to the dataframe - one dataframe!


# APPEND to these 2 dataframes


# Opsim frame, contains OpSim field data (coords, nobs, etc)
OpSimDF = pd.read_csv(path + '0603output_file.csv', header = 0, nrows = 0)#Names for fist lines in datafiles
# Binary dataframe, contains data on binaries in that OpSim field (p, ecc, inc. etc)
BinDF = pd.read_csv(path + '0603output_file.csv', header = 2, nrows = 0)#Names for second data structure in files 




for file in allFiles:
	# Reading in actual data from files, then will append to dataframes created above

	# OpSim data first, should have 1 line of data each
	dat_opsim = pd.read_csv(file, sep = ',', header = 0, nrows = 1)

	# Binaries data read-in - 3rd row on in data files
	dat_bin = pd.read_csv(file, sep = ',', header = 2)
	# print(OpSimDF)
	# Concatenating read-in dataframes to header frames
	OpSimDF = pd.concat(objs = [OpSimDF, dat_opsim])
	BinDF = pd.concat(objs = [BinDF, dat_bin])

# Reindexing OpSim and Binary dataframes, drop = True keeps the old index from being added as a column
OpSimFields = OpSimDF.reset_index(drop = True)
Binaries = BinDF.reset_index(drop = True)

# print(OpSimFields, Binaries)

# Setting up variables from dataframes - pandas Series
PeriodIn = Binaries['p'] # input period -- 'p' in data file - days? - ask Aaron about units
PeriodOut = Binaries['LSM_PERIOD'] #LSM_PERIOD in data file
Mass1 = Binaries['m1']#Binary star #1 mass - Msun
Mass2 = Binaries['m2']#Binary star #2 mass - Msun
radius1 = Binaries['r1']#Binary star #11 radius - Rsun
radius2 = Binaries['r2']#Binary satr #2 radius - Rsun
RA = OpSimFields['OpSimRA']
Dec = OpSimFields['OpSimDec']
OpSimCoords = SkyCoord(RA,Dec,unit=(u.degree, u.degree),frame='icrs')

# N obs in each filter - can use to plot mollweides later 
n_obs_u = OpSimFields['NOpSimObs_u']
n_obs_g = OpSimFields['NOpSimObs_g']
n_obs_r = OpSimFields['NOpSimObs_r']
n_obs_i = OpSimFields['NOpSimObs_i']
n_obs_z = OpSimFields['NOpSimObs_z']
n_obs_y = OpSimFields['NOpSimObs_y']


# ##################################################### Plotting ####################################################################


# Panel plot of mollweides showing Observations per field in each filter 2x3 panel

# Setting up subplots all with mollweide projectsion
f,axarr = plt.subplots(2,3, figsize = (12,6), subplot_kw = dict(projection = 'mollweide'))
f.subplots_adjust(wspace=0.3, hspace = 0)#adjusting subplots

# u filter mollweide
axarr[0,0].grid(True)
scat = axarr[0,0].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
	cmap = 'Blues', c = n_obs_u, vmin = 0, vmax = 200)
axarr[0,0].set_xlabel('u filter')

# g filter mollweide
axarr[0,1].grid(True)
axarr[0,1].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_g, vmin = 0, vmax = 200)
axarr[0,1].set_xlabel('g filter')

# r filter mollweide
axarr[0,2].grid(True)
axarr[0,2].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_r, vmin = 0, vmax = 200)
axarr[0,2].set_xlabel('r filter')

# i filter mollweide
axarr[1,0].grid(True)
axarr[1,0].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_i, vmin = 0, vmax = 200)
axarr[1,0].set_xlabel('i filter')

# z filter mollweide
axarr[1,1].grid(True)
axarr[1,1].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_z, vmin = 0, vmax = 200)
axarr[1,1].set_xlabel('z filter')

# y filter mollweide
axarr[1,2].grid(True)
axarr[1,2].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_y, vmin = 0, vmax = 200)
axarr[1,2].set_xlabel('z filter')


clb = f.colorbar(scat,ax = ax, shrink = 0.75)

plt.show()



















