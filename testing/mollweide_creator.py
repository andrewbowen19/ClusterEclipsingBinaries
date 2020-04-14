# Script to generate Nobs mollweides for colossus and baseline strategies

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
sys.path.insert(0, '/Users/ageller/WORK/LSST/onGitHub/EBLSST/code')

def createMollweide(filename, strategy, create_csv = True):
	'''Function to create mollweide Nobs plots for different OpSim strategies'''
	# Getting OpSim field data and Nobs for each field
	OpS = OpSim()
	OpS.dbFile = f'/Users/andrewbowen/{filename}'
	OpS.getAllOpSimFields()

	# OpSim Coords
	coords = SkyCoord(OpS.RA, OpS.Dec, unit=(u.degree, u.degree), frame='icrs')
	# OpSim field radius
	fieldRad = Angle(1.75, unit = u.degree)

	# Code from Aaron's script
	raGal = coords.icrs.ra.wrap_at(180.*u.degree).degree
	decGal = coords.icrs.dec.wrap_at(180.*u.degree).degree
	lGal = coords.galactic.l.wrap_at(180.*u.degree).degree
	bGal = coords.galactic.b.wrap_at(180.*u.degree).degree

	f, ax = plt.subplots(figsize = (8,5), subplot_kw={'projection': "mollweide"})
	ax.grid(True)
	ax.set_xlabel(r"$RA$",fontsize=16)
	ax.set_ylabel(r"$Dec$",fontsize=16)
	xx = np.where(OpS.Nobs > 0)
	ax.scatter(raGal[xx]*np.pi/180.,decGal[xx]*np.pi/180., s = fieldRad, c=OpS.Nobs[xx], 
	           cmap='tab10', alpha = 1, vmin=0, vmax=1000)
	ax.set_title(strategy, fontsize = 16)
	plt.show()

	if create_csv:
		dat = pd.DataFrame({'OpSim ID' : OpS.ID, 'OpSim RA': OpS.RA, 'OpSim Dec': OpS.Dec, 'Nobs': OpS.Nobs})
		dat.to_csv(strategy + '-OpSim-FieldData.csv')








createMollweide('baseline2018a.db', 'baseline')

