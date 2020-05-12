'''
Script to run analyse for each cluster scenario
Walks through all directories in file tree
Runs analyse script (rewritten as object) for each scenario subdir
Each cluster type dir should have 4 strat/crowding scenario subdirs:
{Baseline-No Crowding, Baseline-Crowding, Col-No Crowding, Col-Crowding}

This should take a while to run
Make sure each scenario sub-directory has updated data files

Pulling LISA candidate WD binaries (P < ~ 2 hours)

Paper here: https://arxiv.org/abs/1907.00014

'''
# Quest version of script


import numpy as np
import pandas as pd
import os
from analyseClusterLISA import analyseClusterLISA
# from analyseCluster import analyse
# from analyseCluster import analyse

# dictionary containig scenario subdir names/acronyms
scenario_dict = {'baseNoCrowd': ['B', 'N'], 'baseCrowd': ['B', 'C'],
                 'colNoCrowd': ['C', 'N'], 'colCrowd': ['C', 'C']}
# Cluster type dictionary
clusterTypes = {'GlobularClusters': 'G', 'OpenClusters': 'O',
                'm10': 'M10', 'm67': 'M67'}

print('Analyzing all clusters including LISA candidate WD binaries...')

# looping through file tree to cover all 16 scenarios
for root, dirs, files in os.walk('./clusters', topdown=True):
    haveStrat = False
    haveCrowd = False
    haveCluster = False
    whiteDwarfs = False
    crowd = None
    stratOpSim = None
    cluster = None

    # Calling analyse for each scenario and cluster type
    for d in dirs:
        # Looping each dir in file tree - calling analyse for each
        path = root + '/' + d  # path to files from this working dir

        # Pulling strat and crowd scenarios from dir names as keys
        if d in scenario_dict:
            print('Have viewing and crowding scenario!', stratOpSim, crowd)
            stratOpSim = scenario_dict[d][0]
            crowd = scenario_dict[d][1]
            haveStrat = True
            haveCrowd = True

        # # Checking cluster type
        if root[11::] in clusterTypes:
            clust_dir = root[11::]
            print('Cluster type acquired!', clusterTypes[clust_dir])
            cluster = clusterTypes[clust_dir]
            haveCluster = True

        # Create cluster analyse object for White Dwarf case
        # Looking for only LISA WDs in this script so no 'searchWD' param
        if haveCrowd and haveStrat:
            print('Running analyse with WDs...')
            ac = analyseClusterLISA(path, cluster, stratOpSim, crowd)
            ac.analyse(path, cluster, stratOpSim, crowd)

    print('')

# TODO: post to quest and run for all clusters
