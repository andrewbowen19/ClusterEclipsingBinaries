""" Using data files from file_compile.py, this should be 
	quicker and easier so we don't have to keep joining dataframes and whatnot 

	Want to now find closest OpSim fields (if there are multiple overlaps)"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle#, search_around_sky
import astropy.units as u

###########################################################################################################################################################

# Reading in gc_data.txt and oc_data.tx files to speed up process, no joining (done in file_compile.py)
names_gc = ['ID_x', 'Name', 'RA', 'DEC', 'L','B','R_Sun','R_gc','X','Y', 'Z', 'key_0','[Fe/H]_x', 'wt', 'E(B-V)_x',\
 'V_HB','(m-M)V_x', 'V_t', 'M_V,t', 'U-B', 'B-V', 'V-R', 'V-I', 'spt', 'ellip', 'ID_y', 'v_r', '+/-', 'v_LSR' ,'sig_v' ,'+/-.1', 'c', 'r_c', 'r_h', 'mu_V',\
  'rho_', 'lg(tc)', 'lg(th)', 'Mcl[Msun]', 'rh[pc]', '[Fe/H]_y', 'age[Gyr]', '(m-M)V_y', 'E(B-V)_y', 'log10(rho[Msun]/pc^3)',\
 'rc', 'sigma0[km/s]', 'esigma0[km/s]', 'fb', 'efb', '[M/H]', 'Rgc[kpc]','Rsun[kpc]']

names_oc = ['name', 'RA', 'Dec', 'l', 'b', 'Dist Mod', 'EB-V', 'Age', 'ST' ,'Z', 'Diam', 'Fe/H', 'MRV',\
 'pm RA', 'pm Dec', 'logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', 'Rhm[pc]',\
  '[Fe/H]K]', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt' ,'Rgc' ,'z' ,'Diam[pc]', 'd[pc]']


path = '/Users/andrewbowen/ceb_project/data/'
GCs = pd.read_csv(path + 'GC_data/gc_data.txt', sep = ' ', header = 0, names = names_gc)
OCs = pd.read_csv(path + 'OC_data/oc_data.txt', sep = ' ', header = 0)
name = OCs['name']

# Dropping NaN values for diameter (for checking for multiple overlaps later)
OCs_nan_drop = OCs.dropna(axis =0, how = 'any', subset = ['name','RA','Dec','Diam'])

OCs_rad = OCs_nan_drop['Diam']/2#radius from WEBDA catalog (most values entered)

# setting up Angle instance for comparing multiple OpSim overlaps later
OCs_rad_angle = Angle(OCs_rad, unit = u.degree)
oc_rc = OCs_rad_angle#consistency with overlap used later (for loops)


# Column name variables for OCs/GCs
oc_diam = OCs['Diam[pc]'] * u.arcmin#OC diameter in parsecs (converted into arcmin)
oc_webda_diam = OCs['Diam'] /2
oc_webda_diam = oc_webda_diam.dropna(how = 'any')

# Definig radii for cluster overlaps later (use of Angle method)
oc_halfdiam = OCs['Diam[pc]'] / 2
oc_halfd = oc_halfdiam.dropna(how = 'any')

oc_rhm = OCs['Rhm[pc]'] *2
oc_rhm = oc_rhm.dropna(how = 'any')
oc_radius = oc_webda_diam.append(oc_rhm)#appending two radii series for open clusters (not every cluster has a listed radius)
oc_radius = oc_radius.dropna(how = 'any')#has all three radii from data table (WEBDA, Piskunov, Solaris)

# Turning into angle radius
oc_rad = Angle(oc_rhm, unit = u.degree)#based off half mass radius (more entries that Diam)
oc_halfdiam = Angle(oc_halfd, unit = u.degree)
# oc_rc = Angle(oc_radius, unit = u.degree)

gc_rhm = GCs['r_h'] * u.arcmin
gc_rad = GCs['rh[pc]'] * 2#using this as the 'true' size of the cluster on the sky (for overlap)
# gc_rad = Angle(gc_rhm, unit = u.arcmin)#creating Angle instance with GC radius

# Dropping duplicate cases in RA/Dec
GCs = GCs.loc[: , ~GCs.columns.duplicated()]
oc_ra = OCs['RA']

OC_RA = oc_ra.dropna(how = 'any')#dropping NaN values in coordinate array for Open Clusters
oc_dec = OCs['Dec']
OC_DEC = oc_dec.dropna(how = 'any')
GC_RA = GCs['RA']#Don't need to drop coords here b/c GCs have coords for all objects (Harris files)
GC_DEC = GCs['DEC']

# Using OpSim data table (RA, Dec, ID for all opsim fields in school), reading in from csv
dat_obs = pd.read_csv('/Users/andrewbowen/ceb_project/data/opsim_field_data.csv', \
		delim_whitespace = True, header =0, names = ['ID', 'RA', 'Dec', 'Nobs'])

# Setting up column names
OpSimID = dat_obs['ID']
OpSimRA = dat_obs['RA']
OpSimDec = dat_obs['Dec']
nobs = dat_obs['Nobs']

# Creating Angle instance within astropy.SkyCoord (turning into degrees)
GC_RA = Angle(GC_RA, unit = u.hour)
OC_RA = Angle(OC_RA, unit = u.hour)
GC_RA = GC_RA.degree
OC_RA = OC_RA.degree

# Creating astropy coord arrays for Open and Globular Clusters
OC_Coords = SkyCoord(OC_RA, OC_DEC, unit=(u.degree, u.degree), frame='icrs')
GC_Coords = SkyCoord(GC_RA, GC_DEC, unit=(u.degree, u.degree), frame='icrs')
lsst_Coords = SkyCoord(OpSimRA, OpSimDec, unit=(u.degree, u.degree), frame='icrs')

# Wrapping coords
lsst_RAWrap = lsst_Coords.ra.wrap_at(180.*u.degree).radian
lsst_DecWrap = lsst_Coords.dec.wrap_at(180.*u.degree).radian

gc_imin, gc_sep2d, gc_dist3d = GC_Coords.match_to_catalog_sky(lsst_Coords)#imin gives array position of matched ID, not id itself
oc_imin, oc_sep2d, oc_dist3d = OC_Coords.match_to_catalog_sky(lsst_Coords)

# Matching coords to OpSim Coords - gives OpSim fields which are closest to clusters
new_gcs = dat_obs.iloc[gc_imin]#using array from gc_imin

# Turning matched coord list into SkyCoords for plotting
gc_match_coords = SkyCoord(new_gcs['RA'], new_gcs['Dec'], unit=(u.degree, u.degree), frame='icrs')

# Open Clusters
OCs['ID'] = pd.Series(oc_imin)
new_ocs = dat_obs.iloc[oc_imin]
oc_match_coords = SkyCoord(new_ocs['RA'], new_ocs['Dec'], unit=(u.degree, u.degree), frame='icrs')

# OpSim Fields: 9.6 deg^2 (circular) - circumference of 3.5 deg (r = 1.75deg)

# Setting up coordinates for multiple check (with missing diams checked)
OC_mult_ra = OCs_nan_drop['RA']
OC_mult_dec = OCs_nan_drop['Dec']
oc_multcheck_coords = SkyCoord(OC_mult_ra, OC_mult_dec, unit=(u.degree, u.degree), frame='icrs')

# Setting up OpSim field radius criterion
r_opsim_arc = Angle((1.76*60), unit = u.arcmin)#width of each opsim field (9.6 deg^2 area, 3.5 deg circumference, 1.75 deg radius) - for GCs, radii given in arcmin
r_opsim = Angle(1.75, unit = u.degree)#for open clusters (in degrees)

# Setting up lists for Globular Cluster for loop
gc_overlap = []
gc_overlap_coord = []
gc_overlap_ra = []
gc_overlap_dec = []
gc_nNear = []

# Globular Cluster For loop - checking overlaps
for index, row in GCs.iterrows():

	coords = SkyCoord(row['RA'], row['DEC'], unit = (u.degree, u.degree), frame = 'icrs')#lsst coords to be compared (all OpSim field coords)
	rad_deg = Angle(row['r_c'], unit = u.arcmin)#cluster radius (half-mass in )
	
	imin, sep1, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 1)#nthneighbor gives multiple OpSim fields
	imin, sep2, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 2)#nthneighbor gives multiple OpSim fields
	imin, sep3, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 3)#nthneighbor to 3rd neighbor

	sep1 = sep1 * 60
	sep2 = sep2 * 60
	sep3 = sep3 * 60

	nNear = 0
	if (sep1-rad_deg < r_opsim_arc):
		nNear += 1
		
	if (sep2-rad_deg < r_opsim_arc):
		nNear += 1
		
	if (sep3-rad_deg < r_opsim_arc):
		nNear += 1

	if (nNear > 1):
		gc_nNear.append(nNear)
		gc_overlap.append(sep1)#could be smarter and give the largest one that overlaps
		gc_overlap_ra.append(row['RA'])
		gc_overlap_dec.append(row['DEC'])
				
	else:
		pass

# Can print numbers of GCs that have overlaps
print("nNear > 1 : ",len(gc_nNear))
print("nNear == 0 : ", len(np.where(nNear == 0)))
print('Total # Globular Clusters:', np.size(GCs['r_c']))

# Globular Cluster coords for multiple overlaps that pass difference criteria
gc_mult_coords = SkyCoord(gc_overlap_ra, gc_overlap_dec, unit = (u.degree, u.degree), frame = 'icrs')

# Setting up lists for open clusters loop
oc_overlap = []
oc_overlap_coord = []
oc_overlap_ra = []
oc_overlap_dec = []
oc_nNear = []

# Now for Open Clusters: need to append coords ----- USING CLUSTER COORDS, NOT OPSIM
for index, row in OCs_nan_drop.iterrows():

	if (row['Diam']):
		coords = SkyCoord(row['RA'], row['Dec'], unit=(u.degree, u.degree), frame='icrs')
		rad = row['Diam']/2#radius from WEBDA catalog (most values entered)

		# setting up Angle instance for comparing multiple OpSim overlaps later
		rc = Angle(rad, unit = u.arcmin)

		imin, sep1, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 1)#nthneighbor gives multiple OpSim fields
		imin, sep2, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 2)#nthneighbor gives multiple OpSim fields
		imin, sep3, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 3)#nthneighbor gives multiple OpSim fields
		#print(ra, dec, sep1, sep2)

		nNear = 0
		if (sep1-rc < r_opsim):
			nNear += 1
		
		if (sep2-rc < r_opsim):
			nNear += 1
		
		if (sep3-rc < r_opsim):
			nNear += 1
		
		# print(row['name'],row['RA'], row['DEC'], rc, sep1, sep2, sep3, sep3-rc, r_opsim, nNear)

		if (nNear > 1):
			oc_nNear.append(nNear)
			oc_overlap.append(sep1) #could be smarter and give the largest one that overlaps
			oc_overlap_ra.append(row['RA'])
			oc_overlap_dec.append(row['Dec'])
				# oc_overlap_coordlist = np.append(oc_overlap_coordlist, coord)#putting coords of separations that match OpSim field criterion
		else:
			pass

# print('Done with Webda loop!!!')

# Dropping rows with no radii - can't compare without radii
OCs_rhm_drop = OCs.dropna(axis =0, how = 'any', subset = ['name','RA','Dec','Rhm[pc]'])
oc_rhm_angle = []

# Converting radii units (from parsecs given by Piskunov)
for r in OCs_rhm_drop['Rhm[pc]']:
	(r * u.parsec).to(u.arcmin, equivalencies = u.parallax())
	oc_rhm_angle.append(r)


pisk_overlap = []
pisk_overlap_ra = []
pisk_overlap_dec = []
pisk_nNear = []

# Piskunov half-mass radius loop:
for index, row in OCs_rhm_drop.iterrows():

	if (row['Rhm[pc]']):
		coords = SkyCoord(row['RA'], row['Dec'], unit=(u.degree, u.degree), frame='icrs')
		rad = row['Rhm[pc]']#radius from WEBDA catalog (most values entered)

		# setting up Angle instance for comparing multiple OpSim overlaps later
		rc = Angle(rad, unit = u.arcmin)
		# print(rc)
		imin, sep1, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 1)#nthneighbor gives multiple OpSim fields
		imin, sep2, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 2)#nthneighbor gives multiple OpSim fields
		imin, sep3, dist3d = coords.match_to_catalog_sky(lsst_Coords, nthneighbor = 3)#nthneighbor gives multiple OpSim fields
		#print(ra, dec, sep1, sep2)

		nNear = 0
		if (sep1-rc < r_opsim):
			nNear += 1
		
		if (sep2-rc < r_opsim):
			nNear += 1
		
		if (sep3-rc < r_opsim):
			nNear += 1
		
		# print(row['name'],row['RA'], row['DEC'], rc, sep1, sep2, sep3, sep3-rc, r_opsim, nNear)

		if (nNear > 1):
			pisk_nNear.append(nNear)
			pisk_overlap.append(sep1) #could be smarter and give the largest one that overlaps
			pisk_overlap_ra.append(row['RA'])
			pisk_overlap_dec.append(row['Dec'])
				# oc_overlap_coordlist = np.append(oc_overlap_coordlist, coord)#putting coords of separations that match OpSim field criterion
		else:
			pass



# Printing statistics for OC multiple overlaps
print(' ')
print('Open Cluster')
print("nNear > 1 : ",len(pisk_overlap))
print("nNear == 0 : ", len(np.where(nNear == 0)))
print("Total OCs with radius:", len(OCs_nan_drop['Diam']))

# Webda multiple overlaps 
oc_overlap_ra = Angle(oc_overlap_ra * 15, unit = u.hour)
oc_overlap_dec = Angle(oc_overlap_dec * 15, unit = u.degree)
oc_multi_coords = SkyCoord(oc_overlap_ra, oc_overlap_dec, unit=(u.degree, u.degree), frame='icrs')

# Setting up angle instances for coordinates later - Piskunov 
piskra = Angle(pisk_overlap_ra, unit = u.degree)
piskdec = Angle(pisk_overlap_dec, unit = u.degree)

# Making instance of OpSim fields with multiple overlapping clusters
pisk_mult_coords = SkyCoord(pisk_overlap_ra,pisk_overlap_dec, unit = (u.degree, u.degree), frame = 'icrs')

############################################################ Plotting ####################################################################################################

# Radius histogram (Globular Clusters)
# plt.hist((gc_rhm.dropna(how = 'any')), bins = 50)#Dropping nan values for radii for GCs
# plt.title('Globular Cluster radii')
# plt.xlim((0,20))
# plt.ylabel('# of binaries')
# plt.xlabel('radii (arcmin)')
# plt.savefig('/Users/andrewbowen/ceb_project/plots/gc_radius_hist.pdf')

# # Open Clusters radii histogram
# plt.hist((oc_radius*u.arcmin), bins = 1500)
# plt.title('Open Cluster radii')
# plt.xlim((0,20))
# plt.ylabel('# of binaries')
# plt.ylabel('radii (arcmin)')
# plt.savefig('/Users/andrewbowen/ceb_project/plots/oc_radius_hist.pdf')
# plt.show()

# Globular Clusters - nNear histogram
# f, ax = plt.subplots(figsize = (8,5))
# ax.hist(gc_nNear, bins = 3, color = 'red')
# ax.set_title('nNear data - GCs')
# ax.set_ylabel('Nbin')
# f.savefig('/Users/andrewbowen/CEB_Project/data/plots/overlaps/gc_nNear_hist.pdf')

# # Open Clusters - WEBDA nNear histogram
# f, ax = plt.subplots(figsize = (8,5))
# ax.hist(oc_nNear, bins = 3, color = 'blue')
# ax.set_title('nNear data - OCs')
# ax.set_ylabel('Nbin')
# f.savefig('/Users/andrewbowen/CEB_Project/data/plots/overlaps/oc_nNear_hist.pdf')

# # Piskunov nNear histogram (OCs)
# f, ax = plt.subplots(figsize = (8,5))
# ax.hist(pisk_nNear, bins = 3, color = 'yellow')
# ax.set_title('nNear data - OCs (webda)')
# ax.set_ylabel('Nbin')
# f.savefig('/Users/andrewbowen/CEB_Project/data/plots/overlaps/webda_nNear_hist.pdf')
# plt.show()

# Cluster radii for plotting
ocDiam = oc_diam.degree
gcDiam = gc_rhm.degree

# Open Cluster Mollweide (with lsst fields)
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap, s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)#Plotting OpSim Fields (all!)
plt.scatter((OC_Coords.ra).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*ocDiam, c = '#A62B1F')
plt.scatter((OC_Coords.ra).wrap_at(180.*u.degree).radian, OC_Coords.dec.radian, s=5*oc_rhm, c = '#A62B1F')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Open Clusters')
plt.savefig('/Users/andrewbowen/ceb_project/plots/oc-coords-poster.pdf', bbox_inches='tight')

# Globular Cluster Mollweide (with lsst fields)
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap, s=5, cmap = 'Blues', c = nobs, vmin = 0, vmax = 250)
plt.scatter((GC_Coords.ra).wrap_at(180.*u.degree).radian, GC_Coords.dec.radian, s=5*gcDiam, c = '#A62B1F')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Globular Clusters')
plt.savefig('/Users/andrewbowen/ceb_project/plots/gc-coords-poster.pdf', bbox_inches='tight')

# GCs match plot
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap,s=5, c = '#5687A6')#colors work for our poster (Dali red-blue)
plt.scatter((gc_match_coords.ra).wrap_at(180.*u.degree).radian, gc_match_coords.dec.radian, s=5, c = '#A62B1F')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Globular Cluster closest OpSim fields')
plt.savefig('/Users/andrewbowen/ceb_project/plots/gc_opsim.pdf', bbox_inches='tight')

# OCs match plot
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap,s=5, c = '#5687A6')
plt.scatter((oc_match_coords.ra).wrap_at(180.*u.degree).radian, oc_match_coords.dec.radian, s=5, c = '#A62B1F')#blue background - red for OCs
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Open Cluster closest OpSim fields')
plt.savefig('/Users/andrewbowen/ceb_project/plots/oc_opsim.pdf', bbox_inches='tight')


# Globular Clusters Multiple
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap,s=5, c = 'blue')
plt.scatter(((gc_mult_coords.ra)*15).wrap_at(180.*u.degree).radian, gc_mult_coords.dec.radian, s=5, c = 'red')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Globular Cluster closest OpSim fields - multiple')
plt.savefig('/Users/andrewbowen/ceb_project/plots/gc-multoverlap.pdf')

# Comparing with these coords
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap,s=5, c = '#13294B')
plt.scatter(((oc_multcheck_coords.ra)*15).wrap_at(180.*u.degree).radian, oc_multcheck_coords.dec.radian, s=5, c = '#7BAFD4')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Open Cluster closest OpSim fields - multiple')


# Plot of piskunov and webda multiple overlaps overlaid - compare with 
plt.figure()
plt.subplot(111, projection = 'mollweide' )
plt.grid(True)
plt.scatter(lsst_RAWrap, lsst_DecWrap, s=5, c = 'blue')
plt.scatter(((oc_multi_coords.ra)).wrap_at(180.*u.degree).radian, oc_multi_coords.dec.radian, s=5, c = 'red')
plt.scatter(((pisk_mult_coords.ra)*15).wrap_at(180.*u.degree).radian, pisk_mult_coords.dec.radian, s=5, c = 'white')#, marker = 'X')
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\delta$')
plt.title('Piskunov & WEBDA overlaps - Open Clusters')
plt.savefig('/Users/andrewbowen/ceb_project/plots/pisk-webda-overlap.pdf')


plt.show()






