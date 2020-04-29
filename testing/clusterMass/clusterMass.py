'''
Script to estimate total cluster mass
Using EBLSST data files compiled by Aaron
'''

import numpy as np
import pandas as pd


gc = pd.read_csv('./data/GCdataForEBLSST.csv')
oc = pd.read_csv('./data/OCdataForEBLSST.csv')

gcMass = gc['mass[Msun]']
ocMass = oc['mass[Msun]']

gcMassSum = np.sum(gcMass)
ocMassSum = np.sum(ocMass)


# Printing results
print(f'The total mass for Globular Clusters is: {gcMassSum} Msol')
print(f'The total mass for Open Clusters is: {ocMassSum} Msol')
print('#######################')
print(f'The total mass for Globular Clusters is: {np.log10(gcMassSum)} log-Msol')
print(f'The total mass for Open Clusters is: {np.log10(ocMassSum)} log-Msol')


