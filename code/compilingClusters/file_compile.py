# Cluster EB Project

# Marin-French: GCs
# Salaris, Piskunov: OCs

"""Parameters we want:
	name, RA, Dec, distance from galactic center, distance from Earth, size (ideally a half-mass radius), 
	number of stars, age, reddening (Av or E(B-V)), metallicity, central density, central velocity dispersion, 
	concentration parameter"""

from pandas import DataFrame
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance, match_coordinates_sky
import glob


# Making .txt file for all data files in directory
path = '/Users/andrewbowen/ceb_project/data/'
filenames = ['Piskunov2008_M.Kharchenko2013_Rhm_FeH.txt',\
	'Solaris2004_viaWEBDA_plusvandenbergh2006_diam_dist.txt', "WEBDA-OC-table-June2013_DavidJames.txt", "mwgc1.txt", "mwgc2.txt", "mwgc3.txt"]


############################# Open Clusters #################################################################################
# WEBDA data file (2013)
webda_df = pd.read_fwf(path + "OC_data/WEBDA-OC-table-June2013_DavidJames.txt", widths = [18,14,15,11,9,8,8,8,9,6,9,9,9,7,7,9], header = 0, \
	names = ['name', 'RA', 'DEC' , 'l', 'b', 'Dist', 'Mod', 'EB-V',\
			'Age','ST', 'Z', 'Diam','Fe/H','MRV','pm RA','pm Dec']) #reading in WEBDA-OC-table-June2013_DavidJames.txt datafile
webda_ID = webda_df['name']
webda_ID = webda_ID.str.replace( ' ','_' )#,regex = True)

webda_df['name'] = webda_ID
webda_df['w flag'] = ['W' for x in range(0,len(webda_df))]

# Converting WEBDA Diameter to pc from arcmin using distance (for later)

webda_diam = webda_df['Diam'].values * u.arcmin
# webda diameter conversion to radians
webda_diam_rad = webda_diam.to(u.radian)

webda_dist = webda_df['Dist'].values * u.pc#webda distance (given in pc)

# Conversion
webda_Diam = np.tan(webda_diam_rad) * webda_dist
W_D = pd.Series(webda_Diam)

# Setting values in WEBDA dataframe to pc diameters (linear size) - needed in cluster_muster.ipynb
webda_df['Diam'] = W_D
# print(webda_df['Diam'].to_string())





# Piskunov (2008)
piskunov_df =  pd.read_csv(path + 'OC_data/PiskunovData_withCoords.csv', sep = ',', header = 0,\
	names = ['name','RA','Dec','logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]','rtK[pc]', 'Rhm[pc]', '[Fe/H]K]', 'distanceK[pc]'])
piskunov_df['p flag'] = ['P' for x in range(0,len(piskunov_df))]

# Solaris (2014)
# Had to add last 2 column names from Van Der Bergh (2006): Diam[pc] - radius in parsecs, d[pc] = distance in parsecs
solaris_df = pd.read_csv(path + 'OC_data/SalarisData_withCoords.csv', sep = ',', header = 0,\
	names = ['name', 'RA', 'Dec','deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt', 'Rgc', 'z','Diam[pc]', 'd[pc]'])
solaris_df['s flag'] =['S' for x in range(0,len(solaris_df))]


# Merging Solaris and Piskunov dataframes
PS = piskunov_df.merge(solaris_df, on = 'name', how = 'outer') #works to merge 2 dfs with overlap -> Solaris/Piskunov
# print(PS)
Pisk_RA = PS['RA_x']
Pisk_Dec = PS['Dec_x']
Sol_RA = PS['RA_y']
Sol_Dec = PS['Dec_y']

PS_RA = Pisk_RA.combine_first(Sol_RA)
PS_Dec = Pisk_Dec.combine_first(Sol_Dec)

# Setting PS column values to new combined RA/Dec series
PS['RA'] = PS_RA
PS['Dec'] = PS_Dec

# # Dropping old labels and adding combined series to df
PS = PS.drop(labels = ['RA_x', 'Dec_x','RA_y', 'Dec_y'], axis =1)
# print(PS)


# Merging Solaris/Piskunov df with with WEDBA df
# OCs = webda_df.join(PS.set_index('name'), on = 'name', how = 'outer') #Should give all open clusters
OCs = webda_df.merge(PS, on = 'name', how = 'outer')

# Merging RA/Dec columns into one column

webda_RA = OCs['RA_x']
webda_Dec = OCs['DEC']

Sol_RA = OCs['RA_y']
Sol_Dec = OCs['Dec']

OC_RA = webda_RA.combine_first(Sol_RA)#combining the two series (webda values taken if both present)
OC_Dec = webda_Dec.combine_first(Sol_Dec)

OCs['RA'] = OC_RA

# print(OCs)
OCs = OCs.drop(labels = ['RA_x', 'DEC','RA_y', 'Dec'], axis = 1)

# Resetting declination column to our good merged values
OCs['Dec'] = OC_Dec

new_OC_cols = ['name',  'RA', 'Dec', 'l', 'b', 'Dist', 'Mod', 'EB-V', 'Age', 'ST', 'Z', 'Diam',\
       'Fe/H', 'MRV', 'pm RA', 'pm Dec', 'w flag', 'logM[Msun]', 'rtP[pc]',\
       'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', 'Rhm[pc]', '[Fe/H]K]',\
       'distanceK[pc]', 'p flag', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't',\
       'sigt', 'logt', 'Rgc', 'z', 'Diam[pc]', 'd[pc]', 's flag']
# Resetting OC column order - Coordinates after name
OCs = OCs[new_OC_cols]
needed_cols = ['name', 'RA', 'Dec', 'w flag', 'p flag', 's flag']

# Dropping Rows where we can't find an RA anywhere
indexNames = OCs[(OCs['RA'] == 'NaN') & (OCs['Dec'] == 'NaN')].index
OCs = OCs.drop(indexNames, axis = 0)

print(OCs[1760:1761])
# Missing RA/Dec values (from WEBDA database)
missingRA = []
missingDec = []



############################### Globular Clusters ##########################################################################

# Reading in mwgc data files separately then create one df with all 3 tables
mwgc1 = pd.read_fwf(path + "GC_data/mwgc1.txt", widths = [12,13,13,14,8,8,6,6,6,6,6], header = 0, \
			names = ['ID', 'Name','RA','DEC','L', 'B', 'R_Sun', 'R_gc',  'X','Y','Z']) #widths are column widths in mwgc1 file, using read_fwf bc of funky 'Name' column in theis file
mwgc2 = pd.read_fwf(path + "GC_data/mwgc2.txt", widths = [12,7,4,5,6,6,7,7,7,6,6,6,5,6], header = 0, \
			names = ['ID','[Fe/H]', 'wt','E(B-V)', 'V_HB','(m-M)V', 'V_t', 'M_V,t','U-B','B-V','V-R', 'V-I', 'spt','ellip'])
mwgc3 = pd.read_fwf(path + "GC_data/mwgc3.txt", widths = [12,8,6,9,7,7,9,7,6,6,7,7,8], header = 0, \
			names = ['ID', 'v_r', '+/-', 'v_LSR', 'sig_v','+/-','c','r_c', 'r_h','mu_V', 'rho_', 'lg(tc)','lg(th)'])

# Harris_Marin-French (2009) - different data than mwgc files (same objects) - Globular Cluster
harris_df = pd.read_table(path + "GC_data/Harris_Marin-Franch_GC_Mcl_rh_FeH_age.txt", delim_whitespace = True, header = 0,\
 names = ['ID','Mcl[Msun]','rh[pc]','[Fe/H]','age[Gyr]','(m-M)V','E(B-V)','log10(rho[Msun]/pc^3)','rc','sigma0[km/s]'\
 	,'esigma0[km/s]','fb','efb','[M/H]','Rgc[kpc]',	'Rsun[kpc]'])

# Harris_Marin-French file and mwgc files have same order for globular clusters

# Resetting ID columns in mwgc dfs
ID1 = mwgc1['ID']
ID2 = mwgc2['ID']
ID3 = mwgc3['ID']
Name1 = mwgc1['Name']


# merging mwgc2 and mwgc3 dfs based on ID cols.
ID1 = ID1.astype(str)
ID2 = ID2.astype(str)
ID3 = ID3.astype(str)

# Let's merge!
mwgc = mwgc2.merge(mwgc3, on = ID3) #Merging mwgc2 and mwgc3 dataframes
MWGC = mwgc1.merge(mwgc, left_index = True, right_index = True) #Merging df1 with merged mwgc df


ID = MWGC['ID'] #Isolates ID column for mwgc files
Glob_Clusters = MWGC.merge(harris_df, left_index = True, right_index = True) #Dataframe for all globular clusters

# Fixing any duplicate columns in GC table - passed eye test of 3 clusters, can check more but should be all set
# Glob_Clusters = Glob_Clusters.drop(['ID_y'])
# Glob_Clusters = Glob_Clusters.drop(['key_0'])

GCs = Glob_Clusters.loc[: , ~Glob_Clusters.columns.duplicated()]


# Writing OC and GC dataframes to a file so we don't have to keep re-reading and joining data tables

# Globular clusters
GCs.to_csv(path + 'GC_data/gc_data.txt', sep = ' ')

# # Open Clusters
OCs.to_csv(path + 'OC_data/oc_data.txt', sep = ' ')


# GC column values (online data source): http://physwww.mcmaster.ca/~harris/mwgc.dat

##################################### OpSim Matching #####################################################################################################

# Globular Clusters
gc_dist = GCs['Rsun[kpc]']
gc_metallicity = GCs['[M/H]']
gc_mh = GCs['[M/H]']
gc_ecc = GCs['ellip']
gc_rhl = GCs['r_h']#glob clusters half-light radius
gc_rhm = GCs['rh[pc]']
gc_age = GCs['age[Gyr]']


# Open Clusters
oc_agep = OCs['log(t[yr])K']#OC age from Piskunov
oc_agew = OCs['Age']#Webda age, may need to check units (OCs)
oc_ages = OCs['logt']#Salaris age (OCs)
oc_rhm = OCs['Rhm[pc]']#OC Half-mass radius age 
oc_diam = OCs['Diam[pc]']#OC diameter in parsecs

oc_diam = oc_diam * u.arcmin

# Mollweide plots - Cluster coords
oc_ra = OCs['RA']
OC_RA = oc_ra.dropna(how = 'any')#dropping NaN values in coordinate array for Open Clusters
oc_dec = OCs['Dec']
OC_DEC = oc_dec.dropna(how = 'any')#Same here, just for Dec instead of RA (OCs)

GC_RA = GCs['RA']#Don't need to drop coords here b/c GCs have coords for all objects
GC_DEC = GCs['DEC']

# Creating astropy coord arrays for Open and Globular Clusters
OC_Coords = SkyCoord(OC_RA, OC_DEC, unit=(u.degree, u.degree), frame='icrs')
GC_Coords = SkyCoord(GC_RA, GC_DEC, unit=(u.degree, u.degree), frame='icrs')

new_path = '/Users/andrewbowen/ceb_project/data/'# /projects/p30137/ageller/testing/EBLSST/fix_omega/output_files/' #Path in quest if necessary
allFiles = glob.glob(path + "/*.csv")

# Reading in 'example' to get header
header_file = '/Users/andrewbowen/CEB_Project/data/opsim_data.txt'

# Reading in data file pulled from quest 
dat_obs = pd.read_table(path + 'opsim_field_data.csv', delim_whitespace = True, header = 0, \
	names = ['ID', 'RA', 'Dec', 'Nobs'])

# print(OCs, GCs, dat_obs)

# Definign OpSim Fields
OpSim_RA = dat_obs['RA'].values#numpy representation of RA/DEC
OpSim_DEC = dat_obs['Dec'].values
OpSim_ID = dat_obs['ID']
# nobs = dat_obs['NOpSimObs_r']


# setting up coords from OpSim data file
lsst_Coords = SkyCoord(OpSim_RA, OpSim_DEC, unit=(u.degree, u.degree), frame='icrs')#OpSim coords from data file on Quest

# Metching GCs and OCs to OpSim fields, outputs 3 different arrays (ID, angular separation, and distance from each other (3d))
gc_imin, gc_sep2d, gc_dist3d = GC_Coords.match_to_catalog_sky(lsst_Coords)
oc_imin, oc_sep2d, oc_dist3d = OC_Coords.match_to_catalog_sky(lsst_Coords)


# Merging with OpSim data to find closest lsst fields -> will plot fields later

# Globular Clusters
GCs['ID'] = pd.Series(gc_imin)
new_gcs = GCs
gc_match = pd.merge(dat_obs, new_gcs, on = 'ID')
gc_match_coords = SkyCoord(gc_match['RA'], gc_match['Dec'], unit=(u.degree, u.degree), frame='icrs')


# Open Clusters
OCs['ID'] = pd.Series(oc_imin)
new_ocs = OCs
oc_match = pd.merge(dat_obs, new_ocs, on = 'ID')
oc_match_coords = SkyCoord(oc_match['RA'], oc_match['Dec'], unit=(u.degree, u.degree), frame='icrs')

##################################### Plotting #####################################################################################################

# Closest lsst field mollweide - Open Clusters
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_Coords.ra.wrap_at(180.*u.degree).radian, lsst_Coords.dec.radian,\
	 s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)
plt.scatter(oc_match_coords.ra.wrap_at(180.*u.degree).radian, oc_match_coords.dec.radian,s=5, c='#cc5500')
plt.title('Closest LSST coords - OCs')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')

# Open Cluster mollweide
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_Coords.ra.wrap_at(180.*u.degree).radian, lsst_Coords.dec.radian,\
	 s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_diam, c = '#cc5500')
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_rhm, c = '#cc5500')
plt.title('Open Clusters')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
# plt.savefig('/Users/andrewbowen/CEB_Project/data/data_files/open_mollweide.pdf')

# Closest lsst field mollweide - Globular Clusters
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_Coords.ra.wrap_at(180.*u.degree).radian, lsst_Coords.dec.radian,\
	 s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)
plt.scatter(gc_match_coords.ra.wrap_at(180.*u.degree).radian, gc_match_coords.dec.radian,s=5, c='#cc3433')
plt.title('Closest LSST coords - GCs')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')

# Globular Cluster Mollweide
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_Coords.ra.wrap_at(180.*u.degree).radian, lsst_Coords.dec.radian,\
	 s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)
plt.scatter((GC_Coords.ra*15).wrap_at(180.*u.degree).radian, GC_Coords.dec.radian, s = 5*gc_rhm, c = '#cc3433')
plt.title('Globular Clusters')
plt.legend(('OpSimFields', 'Globular Clusters'))
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
# plt.savefig('/Users/andrewbowen/CEB_Project/data/data_files/globular_mollweide.pdf')

# plt.show()


# Half-light radius (Globular clusters)
f, ax = plt.subplots(figsize = (8,5))
ax.hist(gc_rhl, bins = 50, range = (0,15), color = 'green')
ax.set_xlabel('Half-light radius (arcmin)')
ax.set_title('Globular Clusters')
# f.savefig('/Users/andrewbowen/CEB_Project/data/data_files/GC_halflight_hist.pdf')

# Half-mass radius (OCs)
f, ax = plt.subplots(figsize = (8,5))
ax.hist(oc_rhm, bins = 100, range = (0,15), color = '#99badd', alpha = 0.75)
ax.hist(gc_rhm, bins = 100, range = (0,15), color = '#000080', alpha = 0.5)
ax.legend(('Open Clusters', 'Globular Clusters'))
ax.set_xlabel('Half-mass radius (pc)')
# f.savefig('/Users/andrewbowen/CEB_Project/data/data_files/halfmass_hist.pdf')

# Cluster Age (Open) - Check for GC ages (would be much older than OCs)
f, ax = plt.subplots(figsize = (8,5))
ax.hist(oc_agep, bins = 50, range = (0,15), color = 'blue')
ax.hist(oc_agew, bins = 50, range = (0,15), color = 'blue')
ax.hist(gc_age, bins = 50, range = (0,15), color = 'orange', alpha = 0.75)
ax.set_xlabel('Cluster Age, [Gyr]')


# Open Clusters mollweide (no OpSim fields)
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_diam, c = '#cc5500')
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_rhm, c = '#cc5500')
plt.title('Open Clusters')


# Globular Clusters mollweide of positions (no OpSim fields)
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter((GC_Coords.ra*15).wrap_at(180.*u.degree).radian, GC_Coords.dec.radian, s = 5*gc_rhm, c = '#cc3433')
plt.title('Globular Clusters')


# All Together Now
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_diam, c = '#cc5500')
plt.scatter((OC_Coords.ra*15).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_rhm, c = '#cc5500')
plt.scatter((GC_Coords.ra*15).wrap_at(180.*u.degree).radian, GC_Coords.dec.radian, s = 5*gc_rhm, c = '#cc3433')
plt.legend(('Open Clusters', 'Globular Clusters'))
plt.title('All Clusters')
# plt.savefig('/Users/andrewbowen/CEB_Project/data/data_files/cluster_molleweide.pdf')








