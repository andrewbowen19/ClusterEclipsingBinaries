# Script to concatenate all the long clusters (and GC/OC output files) -- I'm just annoyed I have to deal w 2 files


import pandas as pd
# import matplotlib.pyplot as plt
import os
import numpy as np


header1 = ['OpSimID','OpSimRA','OpSimDec','NOpSimObs_u','NOpSimObs_g','NOpSimObs_r','NOpSimObs_i','NOpSimObs_z','NOpSimObs_y',
'clusterName','clusterMass','clusterDistance','clusterMetallicity','clusterAge','clusterRhm','clusterVdisp']
header2 = ['p','m1','m2','r1','r2','e','i','d','nobs','Av','[M/H]','appMagMean_r','maxDeltaMag','deltaMag_r'
,'eclipseDepthFrac_r','mag_failure','incl_failure','period_failure','radius_failure','eclipseDepth_failure',
'u_LSS_PERIOD','g_LSS_PERIOD','r_LSS_PERIOD','i_LSS_PERIOD','z_LSS_PERIOD','y_LSS_PERIOD','LSM_PERIOD']

def combineLongDF(path, file1, file2):
	clusterName = file1.replace('_output_file.csv', '')

	# first data file
	clusterDat1 = pd.read_csv(path + '/' + file1, nrows = 1,skiprows = 1 ,  names = header1)
	binDat1 = pd.read_csv(path + '/' + file1, skiprows = 2)

	# print(clusterDat1)
	# print(binDat1)

	clusterDat2 = pd.read_csv(path + '/' + file2, nrows =1,skiprows = 1, names = header1)
	binDat2 = pd.read_csv(path + '/'  + file2, skiprows = 2)

	# print(clusterDat2)
	# print(binDat2)

	binDat = pd.concat([binDat1, binDat2])
	print(binDat)

	# writing dfs to a file
	with open(path + '/new_output/' + clusterName + 'output.csv', 'w') as f:
		clusterDat1.to_csv(f)
	with open(path + '/new_output/' + clusterName + 'output.csv', 'a') as f: # .replace('/input_files', '')
		binDat.to_csv(f)




# combineLongDF('./data/', 'M10__output_file.csv', 'M10_2_output_file.csv')


# Need to write this to go through entire file tree
# making mag plots for long clusters
n = 0
for root, dirs, files in os.walk('.', topdown = True):
	haveName1 = False
	haveName2 = False

	for name in files:
		
		if '__output_file.csv' in name:
			print('Have output file #1...')
			haveName1 = True
			name1 = name
		
		if '2_output_file.csv' in name:
			print('Have output file #2...')
			haveName2 = True	
			name2 = name

		if haveName1 and haveName2:
			print('Making combined files for ', name1, name2, ' in file path :', root)
			combineLongDF(root , name1, name2)
			n += 1
			print(f'New output file # {n} made! ')






