# creates dataframe of population data from PhysiBoSS outputs and exports to csv
# Katie Pletz and John Metzcar 2023-2024. 

# Takes a long time!!!! Took 6.5 hours to run serially on ~250 directories. See "loadDataParallel.py" for a faster version.

import pandas as pd                 # for manipulating data
import pcdl                         # physicell data loader library
# import os, sys, subprocess
import math, os, sys, re
import os.path
import time

################## Pathing Variables ########################

# relative path to input directory
# EDIT TO BE LOCATION OF ORIGINAL MODEL FILES
rel_input_dir = 'leukemia_output/'

# relative path to output directory
# EDIT TO BE LOCATION OF WHERE YOU WANT TO STORE THE NEW MODEL FILES
rel_output_dir = 'dataframes/'

########################### PATHING ###########################
# move to PhysiCell root directory - this assumes this script is one folder below the root
# os.chdir("../")
full_path = os.getcwd()
print(full_path)

# INPUTS
input_dir = os.path.join(full_path, rel_input_dir)
print("Input directory: ")
print(input_dir)

# OUTPUTS
output_dir = os.path.join(full_path, rel_output_dir)
print("Output directory: ")
print(output_dir)

# input("Press Enter to continue...\n Press Ctrl + C to exit...")

# make sure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# list of intervention names
# https://stackoverflow.com/questions/800197/how-to-get-all-of-the-immediate-subdirectories-in-python
interventions = [f.name for f in os.scandir(input_dir) if f.is_dir()] 
interventions = sorted(interventions)

# for each intervention
start = time.time()
for i, intervention in enumerate(interventions):
    start_time_in_loop = time.time()
    # create filepath
    path = input_dir + intervention
    print("Path: "+ path)
    print("Intervention: " + str(i))

    # create a timeseries object
    mcds_ts = pcdl.TimeSeries(path, microenv=False,  settingxml=None, graph=False, verbose=False)

    # create dictionary for this intervention
    # keys = time t; values = # live cells at time t
    data = {'intervention': interventions[i]}
    for mcds in mcds_ts.get_mcds_list():
        df_cell = mcds.get_cell_df()
        live_cells =len(df_cell[(df_cell['dead'] == False)])
        data[mcds.get_time()] = live_cells
    
    if (i == 0): df = pd.DataFrame([data])  # create the dataframe
    else: df.loc[len(df)] = data            # append the dictionary to the dataframe
    df.to_csv(output_dir + 'live_cells.csv', index=False)  # save the dataframe to a csv file (save each time - in case process fails)
    print("total time taken this iteration: ", time.time() - start_time_in_loop)

print("total time taken this loop: ", time.time() - start)

df.to_csv(output_dir + 'live_cells.csv', index=False)  # save the dataframe to a csv file

pd.read_csv(output_dir + 'live_cells.csv')  # read the csv file to check that it saved correctly

#     /N/project/pc_patterns/PhysiLab/scripts

#       WARNING: The script f2py is installed in '/N/u/jpmetzca/BigRed200/.local/bin' which is not on PATH.
#   Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.

#     WARNING: The script natsort is installed in '/N/u/jpmetzca/BigRed200/.local/bin' which is not on PATH.
#   Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.