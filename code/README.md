Welcome to the code directory of our Cluster Eclipsing Binaries repository.

This directory features several sub-directories that contain code for various aspects of our project. These are listed below along with instructions for how to best implement the software. Please note that requisite data files may not be included in the GitHub repository, but are available for download via Box.

* **clusterMass**: This directory is for the purpose of estimateing summed and average masses for both Open and Globular clusters. The data files, which include various cluster parameters, reside in the 'data' subdirectory. In order to produce these numbers and print them to the screen run 'clusterMass' via python command line. 
* **whiteDwarfs**: This directory contains the requisite python scripts for our white dwarf analysis. It should be noted that the white dwarf histogram csv files produced by these versions of the analyceCluster scripts reside in this as well as the **flex** directory. In addition, corner plotting scripts for each observing scenario can be found in this directory
* **mollweides** This directory contains code as well as data files for producing the mollweide projection plots included both in our thesis and poster. The csv files needed are both for individual cluster data as well as *OpSim* field observation data. These are located in the **clusters** directory as well as the **data** directory. The class fieldOverlap can be run to check for the number of clusters with multiple overlaps of *OpSim* fields. This can be run for fields with and without observations. In addition, the mollweide_creator script can be run to produce *OpSim* mollweide plots for several observing strategies. Simply comment out the function calls and run from the command line.




