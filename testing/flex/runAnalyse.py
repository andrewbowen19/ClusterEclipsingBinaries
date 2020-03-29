'''
Script to run analyse for each cluster scenario
Walks through all directories/files in our file tree, runs for each scenario subdirectory
Each of the 4 cluster type directories should have 4 strategy/crowding scenario subdirectories:
{Baseline-No Crowding, Baseline-Crowding, Colossus-No Crowding, Colossus-Crowding}

This should take a while to run (ideally we should copy this directory to Quest and run it)
Make sure each scenario sub-directory has updated data files -- will need Quest paths when uploaded to Quest

'''



import numpy as np
import pandas as pd
import os
import analyseCluster
# from analyseCluster import analyse

# dictionary containig scenario subdir names and they're corresponding inputs as params for analyseCluster
scenario_dict = {'baseNoCrowd' : ['B','N'], 'baseCrowd' : ['B','C'], 'colNoCrowd' : ['C','N'], 'colCrowd' : ['C','C']} 
clusterTypes = {'GlobularClusters':'G', 'OpenClusters':'O', 'm10':'M10', 'm67':'M67'}

# looping through file tree to cover all 16 scenarios
for root, dirs, files in os.walk('./clusters', topdown = True):
	haveStrat = False
	haveCrowd = False
	haveCluster = False
	whiteDwarfs = False

		# Now need to instantiate analyseCluster for each scenario (clusterType, baseline/colossus, crowd/no-crowd)
	for d in dirs:
		# Running through each directory in file tree
		path = root + '/' + d # path to files from this script's home directory
		print('path', path)
		
		# Pulling strategy and crowding scenarios from directory names -- tying to scenario_dict
		if d in scenario_dict:
			print('Have scenario!')
			stratOpSim = scenario_dict[d][0]
			crowd = scenario_dict[d][1]
			haveStrat = True
			haveCrowd = True

		# Checking cluster type
		if d in clusterTypes:
			print('Cluster type acquired!')
			cluster = clusterTypes[d]
			haveCluster = True

		# Checking if white dwarf case is being created
		if 'wd' in (root + d):
			whiteDwarfs = True
		print('')
		
		# Create cluster analyse object without White Dwarf case
		if haveCluster and haveCrowd and haveStrat:
			print('Running analyse without WDs...')
			run = analyseCluster(cluster, stratOpSim, crowd, False)
			run.analyse(cluster, stratOpSim, crowd, False)

		# create analyse object for White Dwarf case
		if haveCluster and haveCrowd and haveStrat and whiteDwarfs:
			print('Running analyse WITH WDs...')
			run = analyseCluster(cluster, stratOpSim, crowd, True)
			run.analyse(cluster, stratOpSim, crowd, True)
	print('')






















