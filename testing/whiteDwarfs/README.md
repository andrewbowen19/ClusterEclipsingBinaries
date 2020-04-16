## WHITE DWARF ECLIPSING BINARIES ##
This subdirectory of our ClusterEclipsing Binaries project is for our analysis of White Dwarf (WD) Binaries in our samples

These were first noticed in some of our population-scenario corner plots
We have since identified some WDs present across various viewing scenarios and recovery sub-populations

Our class and subsequent methods in whiteDwarfBinary selects WD pairs from our analyse script output files.
We will pull these from our output .csv files from our flex directory (produced by lisaRunAnalyse)

Some of the WDs in our sample will be detectable by both LSST and LISA. These are marked in our LISA files ('histDataLISA.csv')
General WD pairs (possibly NOT detectable by LISA but only LSST) are present in this directory and subdirectories

In addition, plots created in this subdirectory are of binaries including 1 or 2 WDs. We primarily concern ourselves with the WD binaries in this analysis.

In our plots folder, the combined plots are divided into plot type: hist and scatter.
These 'combined plots' show scatter plots and hists across all 4 observing scenarios (base-Crowd, base-No Crowd, colossus-No Crowd, colossus-Crowd)

*Note:* the csv files with the extensions 'wdBinaryPairs.csv' are pulled from the 'WD-histData.csv' files in the 
clusters' subdir - those files were produced from the wd-analysescripts located in those subdirectories of clusters

