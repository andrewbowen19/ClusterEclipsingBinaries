# Adding Salaris 2004 coords to the data we have - want RA/Dec for every Salaris Cluster


import pandas as pd

path = '/Users/andrewbowen/ceb_project/data/OC_data/'

# names from new salaris data table
new_names = ['N','Identifier','Otype','RA','Dec','Mag U','Mag B','Mag V','Mag R','Mag I','Sp type','#ref 1850 - 2019','#notes']

# old salaris column names, from: Solaris2004_viaWEBDA_plusvandenbergh2006_diam_dist.txt
sol_names = ['name', 'deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt', 'Rgc', 'z','Diam[pc]', 'd[pc]']

# Reading in 2 Salaris files
new_sol = pd.read_table(path + 'new_salaris_coords.txt', sep = '\t',  header = 0, names = new_names)
old_sol = pd.read_table(path + 'Solaris2004_viaWEBDA_plusvandenbergh2006_diam_dist.txt', delim_whitespace = True, \
	header = 0, names = sol_names)

sol_RA = new_sol['RA']
sol_dec = new_sol['Dec']

# Maybe try to dump excess stuff and merge on names column?

# print(new_sol['RA'])
new_names = ['N','name','Otype','RA','Dec','Mag U','Mag B','Mag V','Mag R','Mag I','Sp type','#ref 1850 - 2019','#notes']
new_sol.columns = new_names
new_sol = new_sol[new_names]#resetting column names


# Making names standrard across both files (from file_compile)

New_Sol_Names = new_sol['name']
New_Sol_Names = New_Sol_Names.str.replace(' ', '_')
new_sol['name'] = New_Sol_Names#putting it back into the df

# merging new and old salaris dfs
all_sol = old_sol.join(new_sol.set_index('name'), on = 'name', how = 'outer')
# print(all_sol.columns)

# List of columns to drop from all_sol
dropped_cols = ['N','Otype','Mag U','Mag B', 'Mag V', 'Mag R', 'Mag I', 'Sp type', '#ref 1850 - 2019','#notes']
Salaris_df = all_sol.drop(labels = dropped_cols, axis = 1)


# Final columns to use for our table
final_salaris_cols = ['name', 'RA', 'Dec','deltaV', 'sigdV', '[FeH]', 'sigFeH', 't', 'sigt', 'logt',\
       'Rgc', 'z', 'Diam[pc]', 'd[pc]']
Salaris_df = Salaris_df[final_salaris_cols]#rearrangin column order

# Sending new df to csv file - will need to fill in missing 10 RA/Dec values manually
Salaris_df.to_csv(path + 'SalarisData_withCoords.csv', sep = ',', header = final_salaris_cols)


