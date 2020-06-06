# Piskunov Coords script (same idea as salaris_coords) - want to add coordinates from Simbad to Piskunov table

import pandas as pd

path = '/Users/andrewbowen/ceb_project/data/OC_data/'

# Reading in old piskunov
piskunov_df =  pd.read_table(path + 'Piskunov2008_M.Kharchenko2013_Rhm_FeH.txt', delim_whitespace = True, header = 0, \
	names = ['name','logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]','rtK[pc]', 'Rhm[pc]', '[Fe/H]K', 'distanceK[pc]'])

# Reading in new Pisk table
new_names = ['N','Identifier','Otype','RA','Dec','Mag U','Mag B','Mag V','Mag R','Mag I','Sp type','#ref 1850 - 2019','#notes']

new_pisk = pd.read_table(path + 'Piskunov_NewCoords.txt', sep = '\t',  header = 0, names = new_names)

# Reordering column names in new df (with coords)
new_names = ['N','name','Otype','RA','Dec','Mag U','Mag B','Mag V','Mag R','Mag I','Sp type','#ref 1850 - 2019','#notes']
new_pisk.columns = new_names
new_pisk = new_pisk[new_names]

# Adding underscore to new names so merge can work
New_Pisk_Names = new_pisk['name']
New_Pisk_Names = New_Pisk_Names.str.replace(' ', '_')
new_pisk['name'] = New_Pisk_Names#putting it back into the df


# Merging dfs
all_pisk = piskunov_df.join(new_pisk.set_index('name'), on = 'name', how = 'outer')
# print(all_sol.columns)

# List of columns to drop from all_sol
dropped_cols = ['N','Otype','Mag U','Mag B', 'Mag V', 'Mag R', 'Mag I', 'Sp type', '#ref 1850 - 2019','#notes']
Piskunov_df = all_pisk.drop(labels = dropped_cols, axis = 1)

final_pisk_cols = ['name', 'RA', 'Dec','logM[Msun]', 'rtP[pc]', 'log(t[yr])K', 'rcK[pc]', 'rtK[pc]', \
	'Rhm[pc]', '[Fe/H]K', 'distanceK[pc]']
Piskunov_df = Piskunov_df[final_pisk_cols]

# print(Piskunov_df)

Piskunov_df.to_csv(path + 'PiskunovData_withCoords.csv', sep = ',', header = final_pisk_cols)