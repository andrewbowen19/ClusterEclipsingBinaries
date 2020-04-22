This testing directory is for us to re-create the Nobs mollweide plots we used in our AAS poster

To fully utilize OpSim, you would have to re-download the baseline and colossus db files
**Note**: These files take up large amounts of space and need the OpSim.py file to be opened

Luckily, we have exported the relevant data from those files (Nobs, RA, Dec, Field ID) into the csv files located in our 'output' subdir
Use these csv files in the future when recreating the mollweides (saved as .pdf and .png files in *plots*)

May want to create a 'plotting' class for this project that we can just call as an import package in other scripts

UPDATE: Was using incorrect colossus download file (colossus_2665 was baseline-style; incorrect)
- Now using colossus 2664, which actually weighs the galactic plane

OpSim File sql db downlaod link: http://astro-lsst-01.astro.washington.edu:8080/?runId=16


