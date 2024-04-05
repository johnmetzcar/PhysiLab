
#!/usr/bin/env python
import pandas as pd                 # for manipulating data
import pcdl                         # physicell data loader library
# import os, sys, subprocess
import math, os, sys, re
import os.path
import time


# MPI stuff
from mpi4py import MPI 
import numpy as np
import os, sys, subprocess

def loadData(intervention, path):
    path = path + intervention
    print("Path: "+ path)
    print("Intervention: " + intervention)

    # create a timeseries object
    mcds_ts = pcdl.TimeSeries(path, microenv=False,  settingxml=None, graph=False, verbose=False)

    # create dictionary for this intervention
    # keys = time t; values = # live cells at time t
    data = {'intervention': intervention}
    for mcds in mcds_ts.get_mcds_list():
        df_cell = mcds.get_cell_df()
        live_cells =len(df_cell[(df_cell['dead'] == False)])
        data[mcds.get_time()] = live_cells

        df = pd.DataFrame([data]) # will end up having to remove the time - but this would allow us to plot each one if we want!!

        df.to_csv(path + 'live_cells.csv', index=False)  # save the dataframe to a csv file


if __name__ == '__main__':
    # PhysiCell_Model, Samples, Replicates = args_run_simulations(sys.argv[1:]) 

    loadData(sys.argv[1], sys.argv[2])
    # # make_output_directories(Samples) --> in create_files.py now

    # NumSimulations = len(Samples) # - which is a list of strings number of samples - this is in fact the number of simulations your run - 
    #         # ID of each aparmaeter set. Would need the right folder for storage and EVERYTHING is ready to run. And the physicell executable is in the right place.
    #         # also need to chagne the seed. 
    #         # I THINK samples is indexless. 
    # NumSimulationsPerRank  = int(NumSimulations/size)

    # # print(NumSimulationsPerRank)
    
    # print(f"Total number of simulations: {NumSimulations} Simulations per rank: {NumSimulationsPerRank}")

    # data = np.array_split(np.arange(NumSimulations),size, axis=0) # split [0,1,...,NumSimulations-1] in size arrays equally (or +1 in some ranks)

    # print(data)