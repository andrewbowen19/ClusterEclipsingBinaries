##Testing directory to estimate summed cluster statistics

Want to include these numbers in thesis
Will use to justify why GCs have higher recovery #s

Pulling data from O/GCdataForEBLSST csv files
*Note:* These use mean masses in pace of NaN values for clusters with no masses listed.

In order to gnerate statistics for different cluster parameters, run *clusterStats* via command line.
The parameters that can be pulled from the cluster data csv files located in the data subriectory are:
* Cluster Age [Myr]
* Cluster velocity dispersion [km/s]
* CLuster distance [pc]
* Cluster half-mass radius [pc]
* Cluster mass [M_sol]
* Cluster metallicity

This script can produce different plot types of cluster parameters for different cluster types as well as different cluster parameters.
The input arguments for the clusterStats object:

* clusterParams: list-like including keys for different cluster parameters (more onkeys later)
* clusterType: string of cluster types; Globular, Open, all (combined across cluster types)

**Keys**:
These shorter strings are used to index into the cluster dataframes
These keys are used when sending parameter arguments to the clusterStats class

key | column label
---- | -------------
'dist' | 'dist[pc]'
'rhm' | 'rhm[pc]'
'mass' | 'mass[Msun]'
'age' | 'age[Myr]'
'z' | '[Fe/H]'
'sigma' | 'sigma_v0_z[km/s]'

So for example, when instantiating the clusterStats object, 'age' would be passed as an argument to grab the 'age[Myr]' column from the cluster dataframe.

It should also be noted that the clusterMass script is deprecated. It's functionality is contained in the much better clusterStats object. Use that instead.
Plots produced from generating this object for each cluster type and parameter can be downloaded [here.](https://northwestern.box.com/s/3a7eb53lepcqtkp1ozrvwkiwpk87binu)
