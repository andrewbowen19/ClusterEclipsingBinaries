Directory for us to make analyse script into a more flexible class.
Want to be able to import it from wherever

Files copied are M10 baseline no-crowding output_files as a test
Making this as a directory with file structure similar to that of WDs

Want this to be able to be called in other scripts
so we can run it in an os.walk loop

can use 'runAnalyse' script to run analyse for each subscenario (16 in total -- OC/GC/LongCluster, Baseline/Colossus, Crowding/No-Crowding)
Need to update directories with proper input files
Will upload folder to Quest and run on there - will make OCs and GCs use less memory


UPDATE: successfully ran lisaRunAnalyse and RunAnalyse script to produce output files for both LISA-candidate WDs as well as circularized binary output files

Output csv files from Quest are posted in this repo. These are for histogram/corner plot data for MS stars as well as WDs (LISA WDs specifically).

There are 2 different histogram csv files produced from our analyseClusterLISA and analyseCluster script
As their names imply, the analyseClusterLISA script has all the functionality of analyseCluster,
while also selecting LISA-LSST WD candidates and producing those histogram files (-histDataLISA.csv files)

In order to use the analyseCluster class for each observing scenario, run either 'runAnalyse' or 'lisaRunAnalyse' depending on if you'd like to include WD analysis

Binary population data for each cluster can be downloaded [here](https://northwestern.box.com/s/wb4cyw9ihne1lffkov4r988l0dk92lvl). The file tree structure of the entire clusters directory is needed for our scripts to run properl. The file tree structure of the entire clusters directory is needed for our scripts to run properly.
