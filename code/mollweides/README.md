## Mollweide plot generation directory

This sub-directory is used to generate mollweide-projection plots used in both [the thesis](https://northwestern.box.com/s/grc92s1vm7qd80jffmnl5guppltlyemt) as well as [the poster presented](https://northwestern.box.com/s/hxuq4fsl76aaxg3m4ctazmjnv6ru56jh) at the AAS conference. 

To generate the mollweide plots from the MySQL files, run the mollweide_creator script once. This will then produce the requisite csv files which contain *OpSim* field coordinates, IDs, and number of observations. These csv files have enough information to produce the same mollweide plots much more quickly. The csv files for different *OpSim* strategies can be downloaded from a box folder [here](https://northwestern.box.com/s/dt40rnlqp33mviszy7mm6b7x87ywoc4r). It is recommended to simply use the csv files to generate the mollweides via the csvMollweide function in mollweide_creator.

OpSim database files for different observing strategies can be downloaded [at this link](http://astro-lsst-01.astro.washington.edu:8080/?runId=16). These files should be located in the __ directory.

* It should be noted that there are variants of different *OpSim* strategies, we utilize the colossus_2664 version, which weighs the galactic plane more heavily. It should be noted that the colossus_2665 is a baseline-style strategy and **should NOT** be used.

* Also note that the nObsMollweides jupyter-notebook is deprectaed, to produce mollweide plots, use the mollweide_creator script

# *OpSim* Field Overlaps

In order to analyse the number of clusters that overlap with multiple *OpSim* fields, run the fieldOverlap python script via command line (<python fieldOverlap.py>). If it is needed to checkfor fields with more than 0 observations, the .loc statement in the variable definitions of <self.fieldRA> and <self.fieldDec>. This script will print to the screen the numbers and percentages of clusters with more than one overlap.

# Overlap Stats

Cluster Type | Percent of Clusters with multiple overlaps
------------ | -----------------------------------------
Globular Clusters | 29.3%
Open Clusters | 33.5%


