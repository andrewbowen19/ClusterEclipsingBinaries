# Script to fix encoding issue with M10 Base-NoCrowd input files

import numpy as np
import pandas as pd


# Only need to read in one header df
header = pd.read_csv('./input_files/M10__output_file.csv', header = 0, nrows = 1)#.drop(columns = 'Unnamed: 0', axis = 1)

# Reading in binary population dfs (both files)
dat1 = pd.read_csv('./input_files/M10__output_file.csv',  skiprows = [0,1])#.drop(columns = 'Unnamed: 0', axis = 1)
dat2 = pd.read_csv('./input_files/M10_2__output_file.csv',  skiprows = [0,1])#.drop(columns = 'Unnamed: 0', axis = 1) #header = 2,


print(dat2.columns)

# Dropping unecessary index column from header and binary dataframes
header = header.drop('Unnamed: 0', axis = 1)
dat1 = dat1.drop('Unnamed: 0', axis = 1)
# dat2 = dat2.drop('Unnamed: 0', axis = 1)


print('writing files...')

# writing df 1 to a file
with open('./input_files/M10__new_output_file.csv', 'w') as f:
	header.to_csv(f, index = False)
with open('./input_files/M10__new_output_file.csv', 'a') as f: # .replace('/input_files', '')
	dat1.to_csv(f, index = False)

# Writing df 2 to a file
with open('./input_files/M10_2__new_output_file.csv', 'w') as f:
	header.to_csv(f, index = False)
with open('./input_files/M10_2__new_output_file.csv', 'a') as f:
	dat2.to_csv(f, index = False)



