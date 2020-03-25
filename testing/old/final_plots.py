# Importing needed libraries
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
# mattplotlib.use('agg')
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance

# Reading in all data files at once
import glob
# path = '/Users/andrewbowen/BO_EB-LSST/code/Opsim'
path ='/projects/p30137/ageller/testing/EBLSST/fix_omega/output_files' # use your path
allFiles = glob.glob(path + "/*.csv")
frame = pd.DataFrame()
frame1 = pd.DataFrame()
frame2 = pd.DataFrame()
list_1 = []
list_2 = []
IDdict = {'RA':np.array([]), 'Dec':np.array([]), 'frac':np.array([])} #For mollewiede plot later

for file_ in allFiles:
    #print(file_)
    dat1 = pd.read_csv(file_, sep = ',', header=0, nrows = 1)
    dat2 = pd.read_csv(file_, sep = ',', header=2)
    dat2['RA'] = dat1['OpSimRA'][0]
    dat2['Dec'] = dat1['OpSimDec'][0]
    list_1.append(dat1)
    #list_2.append(dat2)
    
    if 'p' in dat2.keys():
        #print(dat2['p'].values[0])
        if dat2['p'].values[0] != -1:
            list_2.append(dat2)
            IDdict['RA'] = np.append(IDdict['RA'],dat1['OpSimRA'])
            IDdict['Dec'] = np.append(IDdict['Dec'], dat1['OpSimDec'])
            PeriodIn = dat2['p']
            PeriodOut = dat2['LSM_PERIOD']
            pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values
            f = 0

            recovered_p = dat2['LSM_PERIOD'].loc[pdiff<0.1]
            all_p = dat2['LSM_PERIOD']
#     Putting recovered period into 'frac' lists
            f = recovered_p.shape[0]/all_p.shape[0]


            IDdict['frac'] = np.append(IDdict['frac'],f)
frame1 = pd.concat(list_1) #RA, Dec and id stuff
frame2 = pd.concat(list_2) #Actual data

print('# All Binaries:')
print(frame2.size)

PeriodIn = frame2['p'] # input period -- 'p' in data file
PeriodOut = frame2['LSM_PERIOD'] #LSM_PERIOD in data file
Mass1 = frame2['m1']
Mass2 = frame2['m2']
radius1 = frame2['r1']
radius2 = frame2['r2']
RA = frame1['OpSimRA']
Dec = frame1['OpSimDec']
OpSimCoords = SkyCoord(RA,Dec,unit=(u.degree, u.degree),frame='icrs')

# u filter mollweide
# plt.figure()
f, ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_u'], vmin = 0, vmax = 200)
clbu = plt.colorbar()
f = clbu.set_label(r'Number of Observations', rotation = 270, labelpad = 20)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('u Coordinates', fontsize = 16)
f.savefig('u_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/u_mollweide.pdf')

# g filter mollweide
# plt.figure()
f, ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian,s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_g'], vmin = 0, vmax = 200)
# clb = plt.colorbar()
# clb.set_label(r'Number of Observations', rotation = 270)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('g Coordinates', fontsize = 16)
f.savefig('g_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/g_mollweide.pdf')

# r filter mollweide
# plt.figure()
f, ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_r'], vmin = 0, vmax = 200)
# clb = plt.colorbar()
# clb.set_label(r'Number of Observations', rotation = 270)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('r Coordinates', fontsize = 16)
f.savefig('r_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/r_mollweide.pdf')

# i filter mollweide
# plt.figure()
f, ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s=5,\
            cmap = 'Blues', c = frame1['NOpSimObs_i'], vmin = 0, vmax = 200)
# clb = plt.colorbar()
# clb.set_label(r'Number of Observations', rotation = 270)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('i Coordinates', fontsize = 16)
f.savefig('i_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/i_mollweide.pdf')

# z filter mollweide
# plt.figure()
f, ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian,s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_z'], vmin = 0, vmax = 200)
# clb = plt.colorbar()
# clb.set_label(r'Number of Observations', rotation = 270)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('z Coordinates', fontsize = 16)
f.savefig('z_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/z_mollweide.pdf')

# y filter mollweide
# plt.figure()
f,ax = plt.subplot( projection = 'mollweide' )
ax.grid(True)
ax.scatter(OpSimCoords.ra.wrap_at(180.*u.degree).radian,OpSimCoords.dec.radian, s = 5,\
            cmap = 'Blues', c = frame1['NOpSimObs_y'], vmin = 0, vmax = 200)
# clb = plt.colorbar()
# clb.set_label(r'Number of Observations', rotation = 270)
ax.xlabel(r'$\alpha$', fontsize = 16)
ax.ylabel(r'$\delta$', fontsize = 16)
ax.title('y Coordinates', fontsize = 16)
f.savefig('y_mollweide.pdf') #Save figure
f.savefig('/projects/p30137/abowen/y_mollweide.pdf')


# Period Recovery 
Perr1 = PeriodIn - (0.1*PeriodIn)
Perr2 = PeriodIn + (0.1*PeriodIn)
pdiff = abs(PeriodIn - PeriodOut)/PeriodIn
P_Recover = frame2.loc[pdiff < 0.1]

# ID/All Binaries
Total_Recovered = P_Recover['p'].count
Total_Periods = frame2['p'].count
# print(Total_Recovered)
ID_rate = P_Recover.shape[0]/frame2.shape[0]
print(ID_rate)

# Half-Period correction
Half_Period = 0.5 * PeriodIn
LSM_Half_Period = 0.5 * PeriodOut
hpdiff = abs(Half_Period - LSM_Half_Period)/Half_Period
# HP_Recover = frame2.loc[hpdiff < 0.1]
pdiff = abs(PeriodIn.values - PeriodOut.values)/PeriodIn.values
P_Recover = frame2.loc[(pdiff < 0.1) | (hpdiff < 0.1)]
print(P_Recover)


# Detected Binaries
detected_bins = frame2.loc[PeriodOut != -999]
# print(detected_bins.shape[0])
detection_rate = detected_bins.shape[0]/frame2.shape[0]
print(detection_rate)

# ID/Detected
finder = P_Recover.shape[0]/detected_bins.shape[0]
print(finder)

### Defining Variables for grid plot  ###
# Period
p = frame2['p']
d_p = detected_bins['p']
i_p = P_Recover['p']

# Eccentricity
ecc = frame2['e']
d_ecc = detected_bins['e']
i_ecc = P_Recover['e']

# inclination
inc = frame2['i']
d_inc = detected_bins['i']
i_inc = P_Recover['i']

# Mass
mass1= frame2['m1']
mass2 = frame2['m2']
mass_ratio = mass2/mass1

# Fixed Mass
useM = P_Recover.loc[P_Recover['m1'] != -1]
d_useM = detected_bins.loc[detected_bins['m1'] != -1]
i_mass_ratio = useM['m2']/useM['m1']
d_mass_ratio = d_useM['m2']/d_useM['m1']

# Radius ratio

radius1 = frame2['r1']
radius2 = frame2['r2']
radius_ratio = radius2/radius1
#d_radius = detected_bins['r2']/detected_bins['r1']
useR = P_Recover.loc[P_Recover['r1'] != -1]
i_radius_ratio = useR['r2']/useR['r1']
d_useR = detected_bins.loc[detected_bins['r1'] != -1]
d_radius_ratio = d_useR['r2']/d_useR['r1']

################################################################################################################

IDdict = {'RA':np.array([]), 'Dec':np.array([]), 'frac':np.array([])}


print(IDdict)
PeriodOut = frame2['LSM_PERIOD']
PeriodIn = frame2['p']

# Detected Binaries
detected_bins =frame2.loc[PeriodOut != -999]
# print(detected_bins.shape[0])
detection_rate = detected_bins.shape[0]/frame2.shape[0]
print('detection rate:')
print(detection_rate)

# ID/All Binaries
ID_rate = P_Recover.shape[0]/frame2.shape[0]
print('ID rate:')
print(ID_rate)


# Mollweide plot of percent recovered at each RA/Dec
RA = IDdict['RA']
Dec = IDdict['Dec']
Coords = SkyCoord(RA,Dec,unit=(u.degree, u.degree),frame='icrs')

# plt.figure()
f,ax = plt.subplots( subplot_kw={'projection': "mollweide"} )
ax.grid(True)

scat = ax.scatter(Coords.ra.wrap_at(180.*u.degree).radian,Coords.dec.radian, s = 4, cmap = 'Blues', c = IDdict['frac']*100, vmin = 0, vmax = max(0.5))

#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#ax.title('Period Recovery')
clb = f.colorbar(scat,ax = ax, shrink = 0.75) #color bar
ax.set_label(r'% Periods Recovered', fontsize = 18, rotation = 270, labelpad = 20) #label for color bar
ax.set_xlabel(r'$\alpha$', fontsize = 16)
ax.set_ylabel(r'$\delta$', fontsize = 16) 
f.savefig('mollweide_percent_rec.pdf')
f.savefig('/projects/p30137/abowen/mollweide_percent_rec.pdf') #saving with path on Quest
##################################################################################################################

# Big Grid Plot #

f, axarr = plt.subplots(4,5, figsize = (18,12), sharex = 'col')
f.subplots_adjust(wspace=0.3, hspace = 0)

# Column titles
axarr[0,0].set_title('Period', fontsize = 24)
axarr[0,1].set_title('$M_2$/$M_1$', fontsize = 24)
axarr[0,2].set_title('Eccentricity', fontsize = 24)
axarr[0,3].set_title('$R_2$/$R_1$', fontsize = 24)
axarr[0,4].set_title('Inclination', fontsize = 24)
# Row titles
axarr[0,0].set_ylabel('All Binaries', fontsize = 18)
axarr[1,0].set_ylabel('Observable EBs', fontsize = 18)
axarr[2,0].set_ylabel('Recoverable EBs', fontsize = 18)
axarr[3,0].set_ylabel('Cumulative', fontsize = 18)

# Period
axarr[0,0].hist(np.log10(p), bins = 50, range = (0,10), color = '#13294B')
axarr[1,0].hist(np.log10(d_p), bins = 50, range = (0,10), color = '#356897')
axarr[2,0].hist(np.log10(i_p), bins = 50, range = (0,10), color = '#99badd')
axarr[3,0].hist(np.log10(p), bins = 1000, range = (0,10), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(d_p), bins = 1000, range = (0,10), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].hist(np.log10(i_p), bins = 1000, range = (0,10), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,0].set_xlabel('log(Days)')

# Mass Ratio
axarr[0,1].hist(mass_ratio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,1].hist(d_mass_ratio, bins = 50, range = (0,2), color = '#356897')
axarr[2,1].hist(i_mass_ratio, bins = 50, range = (0,2), color = '#99badd')
axarr[3,1].hist(mass_ratio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(d_mass_ratio, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,1].hist(i_mass_ratio, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

# Eccentricity
axarr[0,2].hist(ecc, bins = 50, range = (0,1), color = '#13294B')
axarr[1,2].hist(d_ecc, bins = 50, range = (0,1), color = '#356897')
axarr[2,2].hist(i_ecc, bins = 50, range = (0,1), color = '#99badd')
axarr[3,2].hist(ecc, bins = 1000, range = (0,1), color = '#13294B',density = True,\
             histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(d_ecc, bins = 1000, range = (0,1), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,2].hist(i_ecc, bins = 1000, range = (0,1), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

# Radius Ratio
axarr[0,3].hist(radius_ratio, bins = 50, range = (0,2), color = '#13294B')
axarr[1,3].hist(d_radius_ratio, bins = 50, range = (0,2), color = '#356897')
axarr[2,3].hist(i_radius_ratio, bins = 50, range = (0,2), color = '#99badd')
axarr[3,3].hist(radius_ratio, bins = 1000, range = (0,2), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(d_radius_ratio, bins = 1000, range = (0,2), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,3].hist(i_radius_ratio, bins = 1000, range = (0,2), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)

#Inclination 
axarr[0,4].hist(inc, bins = 50, range = (0,90), color = '#13294B')
axarr[1,4].hist(d_inc, bins = 50, range = (0,90), color = '#356897')
axarr[2,4].hist(i_inc, bins = 50, range = (0,90), color = '#99badd')
axarr[3,4].hist(inc, bins = 1000, range = (0,90), color = '#13294B',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(d_inc, bins = 1000, range = (0,90), color = '#356897',density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].hist(i_inc, bins = 1000, range = (0,90), color = '#99badd', density = True,\
                histtype = 'step', fill = False,cumulative = True)
axarr[3,4].set_xlabel('Degrees')
plt.show()
f.savefig('Poster_grid_plot.pdf')
f.savefig('/projects/p30137/abowen/Poster_grid_plot.pdf')