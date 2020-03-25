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
#path ='/Users/andrewbowen/ceb_project/data/test_files/' # use your path							
path = '/projects/p30137/ageller/testing/EBLSST/output_files/'
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
OpSimDF = pd.read_csv(path + '1588output_file.csv', header = 0, nrows = 0)#Names for fist lines in datafiles
# Binary dataframe, contains data on binaries in that OpSim field (p, ecc, inc. etc)
BinDF = pd.read_csv(path + '1588output_file.csv', header = 2, nrows = 0)#Names for second data structure in files 




for file in allFiles:
	# Reading in actual data from files, then will append to dataframes created above

	# OpSim data first, should have 1 line of data each
	dat_opsim = pd.read_csv(file, sep = ',', header = 0, nrows = 1)

	# Binaries data read-in - 3rd row on in data files
	dat_bin = pd.read_csv(file, sep = ',', header = 2)
	
	# Not reading in files with all -1 values for binary data (still have OpSim data) (bad files)
	if dat_bin['p'][0] != -1:

		# Concatenating read-in dataframes to header frames
		OpSimDF = pd.concat(objs = [OpSimDF, dat_opsim])
		BinDF = pd.concat(objs = [BinDF, dat_bin])

# Reindexing OpSim and Binary dataframes, drop = True keeps the old index from being added as a column
OpSimFields = OpSimDF.reset_index(drop = True)
Binaries = BinDF.reset_index(drop = True)

print(OpSimFields, Binaries)

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

# Period recovery - check for half-period issue we had last year
Perr1 = PeriodIn - (0.1*PeriodIn)#lower bound on input period error
Perr2 = PeriodIn + (0.1*PeriodIn)#upper bound on input period error
pdiff = abs(PeriodIn - PeriodOut)/PeriodIn#
P_Recover = Binaries.loc[pdiff < 0.05]#Binaries where the difference between input and recovered period is less than 5% (of input period)

#recovery rate of period - 0.07 percent (of 120000 binaries)

# Recovered/All Binaries
Total_Recovered = len(P_Recover['p'])
Total_Periods = len(Binaries['p'])
ID_rate = Total_Recovered/Total_Periods


# Half-Period correction - maybe go over this with Aaron
Half_Period = 0.5 * PeriodIn
LSM_Half_Period = 0.5 * PeriodOut
hpdiff = abs(Half_Period - LSM_Half_Period)/Half_Period
pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values
P_Recover = Binaries.loc[(pdiff < 0.1) | (hpdiff < 0.1)]

# Detected Binaries - eclipse able to be seen by LSST - period not necessarily recovered
detected_bins = Binaries.loc[PeriodOut != -999]#Binaries where a period i 
detected_periods = len(detected_bins['p'])
detection_rate = detected_periods/Total_Periods


# Recovered/Detected
finder = P_Recover.shape[0]/detected_bins.shape[0]


print('% of binaries with period recovered out of ALL binaries:', ID_rate * 100, '%')
print('% of detected binaries out of all binaries:', detection_rate * 100)
print('% Recovered binaries out of detected binaries', finder * 100)


# ##################################################### Dataframes #################################################################
# Setting up separate dataframes for All/Observable/Recovered

allBins = Binaries#All binaries
obsBins = detected_bins#observable binaires
recBins = P_Recover#Recovered binaries

# #####################################################Plotting variables ############################################################

# All Binaries first
all_p = allBins['p']#All period
all_m1 = allBins['m1']#All m1
all_m2 = allBins['m2']
all_r1 = allBins['r1']#All radii
all_r2 = allBins['r2']#All radii
all_e = allBins['e']#All eccentricity
all_i = allBins['i']#All inclination
allMassRatio = all_m1/all_m2
allRadRatio = all_r1/all_r2

# Observable Binaries
obs_p = obsBins['p']#obs period
obs_m1 = obsBins['m1']#obs m1
obs_m2 = obsBins['m2']
obs_r1 = obsBins['r1']#obs radii
obs_r2 = obsBins['r2']#obsradii
obs_e = obsBins['e']#obs eccentricity
obs_i = obsBins['i']#obs inclination
obsMassRatio = obs_m1/obs_m2
obsRadRatio = obs_r1/obs_r2

# Recovered Binaries
rec_p = recBins['p']#obs period
rec_m1 = recBins['m1']#obs m1
rec_m2 = recBins['m2']
rec_r1 = recBins['r1']#obs radii
rec_r2 = recBins['r2']#obsradii
rec_e = recBins['e']#obs eccentricity
rec_i = recBins['i']#obs inclination
recMassRatio = rec_m1/rec_m2
recRadRatio = rec_r1/rec_r2

# ##################################################### Plotting ####################################################################


# Panel plot of mollweides showing Observations per field in each filter 2x3 panel

# Setting up subplots all with mollweide projectsion
f,axarr = plt.subplots(2,3, figsize = (12,6), subplot_kw = dict(projection = 'mollweide'))
f.subplots_adjust(right = 0.9,wspace=0.3, hspace = 0)#adjusting subplots

# Making colorbar
cbar_ax = f.add_axes([12,0,0.1,6])
import matplotlib.cm as cm


Normalize = matplotlib.colors.Normalize(vmin=0, vmax=200, clip=False)
mappable = cm.ScalarMappable(norm = Normalize, cmap = 'Blues')
print(mappable)
clb = f.colorbar(mappable, ax = cbar_ax)
f = clb.set_label(r'Number of Observations', rotation = 270, labelpad = 20)

# u filter mollweide
axarr[0,0].grid(True)
axarr[0,0].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
	cmap = 'Blues', c = n_obs_u, vmin = 0, vmax = 200)
axarr[0,0].set_title('u')

# g filter mollweide
axarr[0,1].grid(True)
axarr[0,1].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_g, vmin = 0, vmax = 200)
axarr[0,1].set_title('g')

# r filter mollweide
axarr[0,2].grid(True)
axarr[0,2].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_r, vmin = 0, vmax = 200)
axarr[0,2].set_title('r')

# i filter mollweide
axarr[1,0].grid(True)
axarr[1,0].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_i, vmin = 0, vmax = 200)
axarr[1,0].set_title('i')

# z filter mollweide
axarr[1,1].grid(True)
axarr[1,1].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_z, vmin = 0, vmax = 200)
axarr[1,1].set_title('z')

# y filter mollweide
axarr[1,2].grid(True)
axarr[1,2].scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
    cmap = 'Blues', c = n_obs_y, vmin = 0, vmax = 200)
axarr[1,2].set_title('y')

f.savefig('/projects/p30137/abowen/CEB/plots/old_plots/filters_mollweide.pdf')
# plt.show()


# ########################################### Panel Hist Plot from Last Summer #######################################
# Big Grid Plot - final_plot.py # ##### May need to define some of the variables/change names

f, axarr = plt.subplots(4,5, figsize = (12,8), sharex = 'col')
f.subplots_adjust(wspace=0.3, hspace = 0)


# Column titles
axarr[0,0].set_title('Period', fontsize = 24)
axarr[0,1].set_title('$M_2$/$M_1$', fontsize = 24)
axarr[0,2].set_title('Eccentricity', fontsize = 24)
axarr[0,3].set_title('$R_2$/$R_1$', fontsize = 24)
axarr[0,4].set_title('Inclination', fontsize = 24)
# Row titles
axarr[0,0].set_ylabel('All', fontsize = 18)
axarr[1,0].set_ylabel('Observable', fontsize = 18)
axarr[2,0].set_ylabel('Recoverable', fontsize = 18)
axarr[3,0].set_ylabel('Cumulative', fontsize = 18)

# Period
axarr[0,0].hist(np.log10(rec_p), bins = 50, range = (0,10), color = '#13294B')
axarr[1,0].hist(np.log10(obs_p), bins = 50, range = (0,10), color = '#356897')
axarr[2,0].hist(np.log10(all_p), bins = 50, range = (0,10), color = '#99badd')
axarr[3,0].hist(np.log10(rec_p), bins = 1000, range = (0,10), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(obs_p), bins = 1000, range = (0,10), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(all_p), bins = 1000, range = (0,10), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].set_xlabel('log(Days)')

# Mass Ratio
axarr[0,1].hist(recMassRatio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,1].hist(obsMassRatio, bins = 50, range = (0,2), color = '#356897')
axarr[2,1].hist(allMassRatio, bins = 50, range = (0,2), color = '#99badd')
axarr[3,1].hist(recMassRatio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(obsMassRatio, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(allMassRatio, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

# Eccentricity
axarr[0,2].hist(rec_e, bins = 50, range = (0,1), color = '#13294B')
axarr[1,2].hist(obs_e, bins = 50, range = (0,1), color = '#356897')
axarr[2,2].hist(all_e, bins = 50, range = (0,1), color = '#99badd')
axarr[3,2].hist(rec_e, bins = 1000, range = (0,1), color = '#13294B',density = True,\
             histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(obs_e, bins = 1000, range = (0,1), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(all_e, bins = 1000, range = (0,1), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

# Radius Ratio
axarr[0,3].hist(recRadRatio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,3].hist(obsRadRatio, bins = 50, range = (0,2), color = '#356897')
axarr[2,3].hist(allRadRatio, bins = 50, range = (0,2), color = '#99badd')
axarr[3,3].hist(recRadRatio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(obsRadRatio, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(allRadRatio, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

#Inclination 
axarr[0,4].hist(rec_i, bins = 50, range = (0,90), color = '#13294B')
axarr[1,4].hist(obs_i, bins = 50, range = (0,90), color = '#356897')
axarr[2,4].hist(all_i, bins = 50, range = (0,90), color = '#99badd')
axarr[3,4].hist(rec_i, bins = 1000, range = (0,90), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(obs_i, bins = 1000, range = (0,90), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(all_i, bins = 1000, range = (0,90), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].set_xlabel('Degrees')
#plt.show()

f.savefig('/projects/p30137/abowen/CEB/plots/old_plots/panel_stats.pdf')










