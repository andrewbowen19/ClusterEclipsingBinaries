### Flex code directory

This sub-directory contains the code which analyzes the results of our cluster-binary simulations. These simulation result files reside in the 'input_files'folders in each cluster-scenario sub-directory (e.g. OpenClusters/colCrowd, which corresponds to the Open Cluster-colossus-with crowding observing scenario). 

In order to create PDF/CDF, mollweide recovery plots, and the short-period binary selection histogram csv files, run either the analyseCluster or analyseClusterLISA python scripts via command line: <python analyseCluster.py> or <python analyseClusterLISA.py>. In order to generate these plots and output files for all observing scenarios, utilize the *runAnalyse scripts* via command line, either: <python runAnalyse.py> to create analyseCluster objects or <python lisaRunAnalyse.py> to create analyseClusterLISA objects for each observing scenario. 

Histogram data output csv files from the analyse scripts are located in the 'data' subdirectories for each observing scenario (e.g. 'GlubularClusters/baseCrowd/data/'). Corresponding white dwarf binary output files reside in the 'wd' sub-directories of the 'data' folders. The Box download link provided below contains the needed files structure, which should be placed in a 'clusters' directory in the same directory as the analyse scripts.

There are 2 different histogram csv files produced from our analyseClusterLISA and analyseCluster script
As their names imply, the analyseClusterLISA script has all the functionality of analyseCluster,
while also selecting LISA-LSST WD candidates and producing those histogram files (-histDataLISA.csv files)

In addition, the analyseCluster and analyseClusterLISA scripts produce txt files with recovery statistics for each observing scenario. These are written to the output/stats/ subdirectory.

Binary population data for each cluster can be downloaded [here](https://northwestern.box.com/s/wb4cyw9ihne1lffkov4r988l0dk92lvl). The file tree structure of the entire clusters directory is needed for our scripts to run properly. The file tree structure of the entire clusters directory is needed for our scripts to run properly.

In addition, the corner_plotter script can be run via <pyhton corner_plotter.py> in order to create and save corner plot pdf files for each observing scenario and binary sub-population (*all*, *observable*, or *recovered*). This script loops through each observing scenario sub-directory and uses the histogram csv files in the 'data' sub-directories to generate corner plots for each. The corner plots are then saved to the 'plots/corner_plots' folder located in the working directory.

## Recovery Statistics
# Globular Clusters
Observing Scenario | N_rec | N_obs | % Recovered
------------------ | ----- | ----- | -----------
*baseline* | 687 | 1890 | **36.3%**
*colossus* | 840 | 2377 | **35.3%**
Difference | 153 | 487  | **-1.0%**

# Open Clusters
Observing Scenario | N_rec | N_obs | % Recovered
------------------ | ----- | ----- | -----------
*baseline* | 49 | 148 | **33.1%**
*colossus* | 111 | 286 | **38.8%**
Difference | 62 | 138  | **5.7%**

Overall, *colossus* will recover the periods of ~13% more eclipsing binaries over the 10-year runtime of the Vera Rubin Observatory.



