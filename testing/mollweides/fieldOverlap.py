# New script to check OpSim overlap for clusters
# Want the percentage of clusters that overlap with more than 1 OpSim field
# OCs should overlap more given that they are closer

# We could compare more strategies (I think there is a 'kraken strategy')
# OpSim strategy db file download link
# http://astro-lsst-01.astro.washington.edu:8080/?runId=16

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os 
from astropy import units as u
from astropy.coordinates import Angle, match_coordinates_sky, SkyCoord
# from .mollweide_creator import plotMollweide

class fieldOverlap(object):

    def __init__(self, strategy):
        self.strategy = None

    def plotMollweide(self, x, y, fieldSize, title, colorbar=True):
        '''
        Function to generate mollweide projection plot
        cVal - values for colormap
        from our mollweide_creator code
        '''
        f, ax = plt.subplots(figsize=(8, 5),
                             subplot_kw={'projection': "mollweide"})
        ax.grid(True)
        ax.set_xlabel(r"$RA$", fontsize=16)
        ax.set_ylabel(r"$Dec$", fontsize=16)
        im = ax.scatter(x, y, s=fieldSize,
                        cmap='viridis', alpha=1, vmin=0, vmax=1000)

        ax.set_title(title, fontsize=24)

        if colorbar:
            cbar = f.colorbar(im)
            cbar.ax.set_ylabel(r'Nobs', fontsize=16, rotation=270, labelpad=8)

        f.savefig(f'./plots/{title}-MultipleOverlaps-mollweide.pdf', format='pdf', bbox_inches='tight')
        # plt.show()

    def wrapCoords(self, coords):
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

    def getClusterCoords(self):
        '''Function to get cluster coordinates and convert from hours to degrees (when necessary)'''
        # Reading in cluster data: be sure to use RA & Dec NOT OpSim RA and OpSim Dec
        self.gcData = pd.read_csv(os.path.join('./clusters', 'GCdataForEBLSST.csv'))
        self.ocData = pd.read_csv(os.path.join('./clusters', 'OCdataForEBLSST.csv'))

        # Need cluster radius for separation checks later
        # Converting to angular size
        self.gcRhm = np.ones_like(self.gcData['dist[pc]'].values) * 2. * u.pc# self.gcData['rhm[pc]'].values * u.pc
        self.ocRhm = np.ones_like(self.ocData['dist[pc]'].values) * 2. * u.pc# self.ocData['rhm[pc]'].values * u.pc

        self.gcDist = self.gcData['dist[pc]'].values * u.pc
        self.ocDist = self.ocData['dist[pc]'].values * u.pc

        # Converitng from radians to deg
        # Angular size of clusters (from distance/2pc rhm value calculation): output in radians
        self.gcSize = np.arctan(self.gcRhm/self.gcDist)# * u.rad
        self.ocSize = np.arctan(self.ocRhm/self.ocDist)# * u.rad
        print("Cluster Sizes (rad): ", self.gcSize, '\n',self.ocSize )

        # Converting to degrees (via parallax)
        self.gcRad_deg = (self.gcSize).to(u.deg)#, equivalencies=u.parallax())
        self.ocRad_deg = (self.ocSize).to(u.deg)#, equivalencies=u.parallax())
        print("Cluster Sizes (deg): ", self.gcRad_deg, '\n', self.ocRad_deg)

        # print(self.gcData['RA[hr]'])
        # Converting cluster RA values to degrees from hours
        self.gcRA_deg = Angle(self.gcData['RA[hr]'], unit=u.hour).degree
        self.ocRA_deg = Angle(self.ocData['RA[hr]'], unit=u.hour).degree
        # print('Right Ascensions for clusters: ', self.gcRA_deg, '\n', self.ocRA_deg)

        # setting up declination values
        self.gcDec = Angle(self.gcData['Dec[deg]'].values, unit=u.deg)
        self.ocDec = Angle(self.ocData['Dec[deg]'].values, unit=u.deg)
        # print('Declinations: ', self.gcDec, '\n' ,self.ocDec)

        # Generating SkyCoord option for matching function later
        self.gcCoords = SkyCoord(self.gcRA_deg, self.gcDec, unit=(u.degree, u.degree), frame='icrs')
        self.ocCoords = SkyCoord(self.ocRA_deg, self.ocDec, unit=(u.degree, u.degree), frame='icrs')
        # print('Coordinates for matching: ', self.gcCoords, self.ocCoords)

    def findNeighbors(self, strategy):
        '''
        Function to check for nearest on-sky neighbors for each cluster
        Can run for multiple strategies
        OpSim strategy fields are given in degrees, clusters are listed in hours (conversion time!)
        '''
        # Setting up field coordinates for matching function
        self.fieldDat = pd.read_csv(os.path.join('output/', strategy + '-OpSim-FieldData.csv'))
        self.fieldID = self.fieldDat['OpSim ID']
        self.fieldRA = Angle(self.fieldDat['OpSim RA'].loc[np.where(self.fieldDat['Nobs'] > 0.)], unit=u.degree) # .loc[np.where(self.fieldDat['Nobs'] > 0.)]
        self.fieldDec = Angle(self.fieldDat['OpSim Dec'].loc[np.where(self.fieldDat['Nobs'] > 0.)], unit=u.degree) # .loc[np.where(self.fieldDat['Nobs'] > 0.)]

        self.fieldCoords = SkyCoord(self.fieldRA, self.fieldDec, unit=(u.degree, u.degree), frame='icrs')
        self.fieldRad = Angle(1.75, unit=u.deg)

        # Lists that will contain # of OpSim neighbors for each cluster (OC/GC)
        nOpSimOverlapsGC = np.zeros_like(self.gcData['rhm[pc]']) 
        nOpSimOverlapsOC = np.zeros_like(self.ocData['rhm[pc]'])


        # Searching for closest neighbor OpSim field
        # Will loop through as many neighbors as we want to check and will test if 
        for i in range(1,5):
            print(f'Checking for {i}th closest OpSim field')
            numNeighborsGC = 0
            numNeighborsOC = 0
            # Checking for 
            gc_idx, gc_sep2d, gc_dist3d = match_coordinates_sky(self.gcCoords, self.fieldCoords, nthneighbor=i)
            oc_idx, oc_sep2d, oc_dist3d = match_coordinates_sky(self.ocCoords, self.fieldCoords, nthneighbor=i)

            # 
            # print("Indices of matched coords: ", len(gc_idx), '\n', oc_idx)
            # print("Angular Separations: ", gc_sep2d, '\n', oc_sep2d)
            
            d_gc = gc_sep2d - self.gcRad_deg
            d_oc = oc_sep2d - self.ocRad_deg
            print('Checking if cluster and field overlap...')
            # Checking if separation is wihtin constraints
            for x in range(0,len(d_gc)-1):
                
                if d_gc[x] < self.fieldRad:
                    nOpSimOverlapsGC[x] += 1

            for x in range(0,len(d_oc)-1):
                
                if d_oc[x] < self.fieldRad:
                    nOpSimOverlapsOC[x] += 1
            print('')
        print('# of overlaps for each globular cluster', nOpSimOverlapsGC)
        print('# of overlaps for each open cluster', nOpSimOverlapsOC)
        # Printing clusters with most overlaps (maxes out at 3)
        print('Biggest overlaps (GCs, OCs): ', np.max(nOpSimOverlapsGC),':', 
                self.gcData['ID'][np.argmax(nOpSimOverlapsGC)] , 
                np.max(nOpSimOverlapsOC),':', self.ocData['ID'][np.argmax(nOpSimOverlapsOC)])
        # Checking for multiple fields

        self.gcData['N OpSim Field Neighbors'] = nOpSimOverlapsGC
        self.ocData['N OpSim Field Neighbors'] = nOpSimOverlapsOC

        # Checking percentages of clusters with Nfields > 1
        multipleOverlapPecentGC = len(nOpSimOverlapsGC[np.where(nOpSimOverlapsGC > 1.)])/len(nOpSimOverlapsGC)
        multipleOverlapPecentOC = len(nOpSimOverlapsOC[np.where(nOpSimOverlapsOC > 1.)])/len(nOpSimOverlapsOC)

        # Printing results
        print('########################################')
        print(f'% of GCs that overlap with multiple OpSim fields for {strategy}: ', multipleOverlapPecentGC * 100,
                len(nOpSimOverlapsGC[np.where(nOpSimOverlapsGC > 1.)]),
                len(nOpSimOverlapsGC))
        print(f'% of OCs that overlap with multiple OpSim fields for {strategy}: ', multipleOverlapPecentOC * 100,
                len(nOpSimOverlapsOC[np.where(nOpSimOverlapsOC > 1.)]),
                len(nOpSimOverlapsOC))

        # Plotting clusters with multipl overlaps
        gcMultipleOverlaps = self.gcCoords[np.where(nOpSimOverlapsGC>1.)]
        ocMultipleOverlaps = self.ocCoords[np.where(nOpSimOverlapsOC>1.)]

        # Wrapping coords for plotting
        raGalGC, decGalGC, lGalGC, bGalGC=self.wrapCoords(gcMultipleOverlaps)
        raGalOC, decGalOC, lGalOC, bGalOC=self.wrapCoords(ocMultipleOverlaps)
        
        # GC plots
        self.plotMollweide(raGalGC*np.pi/180., decGalGC*np.pi/180.,
                            self.fieldRad, 'Globular Clusters' + "-" + strategy, False)
        # OC plots
        self.plotMollweide(raGalOC*np.pi/180., decGalOC*np.pi/180.,
                            self.fieldRad, 'Open Clusters' + "-" +  strategy, False)


# Calling object as test
# fo = fieldOverlap('colossus')
# fo.getClusterCoords()
# fo.findNeighbors('colossus')

# Performing same analysis for both strategies
# Don't think this will change much (except maybe remove/add a few with different Nobs)
for s in os.listdir('./output'):
    file = s.replace('-OpSim-FieldData.csv','')
    print(f'Finding overlaps for {file}')
    fo = fieldOverlap(file)
    fo.getClusterCoords()
    fo.findNeighbors(file)
    print('#############################################')
    print('#############################################')







