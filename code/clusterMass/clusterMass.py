'''
Script to estimate total cluster mass
Using EBLSST data files compiled by Aaron
Will output summed cluster masses as well as mean cluster masses
Written as a check for Andrew's thesis
'''

import numpy as np
import pandas as pd
import os

# Reading in cluster data csv files
# need to be located in data directory in same directory as this script
gc = pd.read_csv(os.path.join('./data/', 'GCdataForEBLSST.csv'))
oc = pd.read_csv(os.path.join('./data/', 'OCdataForEBLSST.csv'))

print(os.path.join('./data/', 'GCdataForEBLSST.csv'))

gcMass = gc['mass[Msun]']
ocMass = oc['mass[Msun]']

gcMassSum = np.sum(gcMass)
ocMassSum = np.sum(ocMass)

gcMassMean = np.mean(gcMass)
ocMassMean = np.mean(ocMass)

# Printing results
print(f'The total mass for Globular Clusters is: {gcMassSum} Msol')
print(f'The total mass for Open Clusters is: {ocMassSum} Msol')
print('#######################')
print(f'The total mass for Globular Clusters is: {np.log10(gcMassSum)} log-Msol')
print(f'The total mass for Open Clusters is: {np.log10(ocMassSum)} log-Msol')

# Printing mean masses as well
print('#######################')
print(f'The mean mass for Globular Clusters is: {gcMassMean} Msol', ' | ', "{:e}".format(gcMassMean))
print(f'The mean mass for Open Cluster is: {ocMassMean} Msol', ' | ' , "{:e}".format(ocMassMean))





