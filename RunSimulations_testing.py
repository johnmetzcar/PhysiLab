from mpi4py import MPI 
import numpy as np
import os, sys, subprocess
import math, os, sys, re
# from HPC_exploration import model, args_run_simulations

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Overall quick notes - 

# - Heber suggested using only this file to run the simulations.
# - We need to make the XML completely - this means changing output directors, diffusion constants, and random number seeds
# - WE need to make the directories (see below)
# - We might consider just using this as inspiration and writing our own code.
# - can test on the SICE systems - and just ee what happens. 
# - adapt to use only this one file

# 02.09.2024 - we will just pretend we have all the XMLs. 

# What do I need then?

def generate_list_of_config_files(path_to_config_files: str):

    files = os.listdir(path_to_config_files)
    # print(files)
    # print(path_to_config_files)

    list_of_PC_configs = []

    #### examine all file names in directory and add ones, via string matching, as needed to list of names of files of interest
    for i in range(len(files)):
        if not re.search('\.xml', files[i]):
            continue

        # a NumPy array could be used here. be faster, but I don't
        # expect huge file lists. So I will just sort as I know how to do that ...

        relative_path_config_file = path_to_config_files + '/' + files[i]
        # will work on full pathing later
        # print(file_full_path)
        print(os.path.dirname('PhysiCellModel.py'))

        list_of_PC_configs.append(relative_path_config_file)
    #### Sort file name list - not needed for this application
    # list_of_PC_configs.sort()

    print(list_of_PC_configs)

    return list_of_PC_configs



# make the directories. Could be used to modify XML files as well - if you were feeding in a bunch of bases models. 
# def get_configFilePath(self,sampleID, replicateID):
#     if (self.configFile_folder): 
#         os.makedirs(os.path.dirname(self.configFile_folder), exist_ok=True)
#         return self.configFile_folder+self.configFile_name%(sampleID,replicateID)
#     else:
#         folder = self.get_outputPath(sampleID, replicateID)
#         return folder+self.configFile_name%(sampleID,replicateID)

# def make_output_directories(list_of_PC_configs): --> in create files now
#     for file_name in list_of_PC_configs:
#         # note - will work to get the aboslute path at some point
#         # os.makedirs(os.path.dirname(file_name), exist_ok=True)

#         subdirectory_name = file_name.split('.')[0] # make sure you don't get fanacy with putting `.` in the file name OR path!!!!!!
#         # print(subdirectory_name)
#         subdirectory_name = subdirectory_name.split('/')[1] # will have to be updated for full path!!!
#         # print(subdirectory_name)        
#         os.makedirs('output_HTC/' + subdirectory_name, exist_ok=True)

# Define the PhysiCell execution            
def model(ConfigFile, Executable):
    # Write input for simulation & execute
    callingModel = [Executable, ConfigFile]
    cache = subprocess.run( callingModel,universal_newlines=True, capture_output=True) # stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if ( cache.returncode != 0):
        print(f"Error: model output error! Executable: {Executable} ConfigFile {ConfigFile}. returned: \n{str(cache.returncode)}")
        print(cache.stdout[-200])

### This is what Heber is using now - and of note, if I am going to use it and then remove data, I need to make sure that the data is stored
# somewhere else. And it needs a unique name - thats whats going here. And yes Co-pilote I will sue teh ID and what not to name it. 
# SToring a little bit differently - but still row by row .... Something to think about. 

# def summ_func(OutputFolder,SummaryFile, dic_params, SampleID, ReplicateID):
#     mcds = pcdl.TimeStep('final.xml',OutputFolder, microenv=False, graph=False, settingxml=None, verbose=False)
#     df_cell = mcds.get_cell_df() 
#     live_cells = len(df_cell[ (df_cell['dead'] == False) ] )
#     dead_cells = len(df_cell[ (df_cell['dead'] == True) ] )
#     data = {'time': mcds.get_time(), 'sampleID': SampleID, 'replicateID': ReplicateID, 'live_cells': live_cells, 'dead_cells': dead_cells, 'run_time_sec': mcds.get_runtime()}
#     data_conc = {**data,**dic_params} # concatenate output data and parameters
#     df = pd.DataFrame([data_conc])
#     # remove replicate output folder
#     rmtree( OutputFolder )
#     df.to_csv(SummaryFile, sep='\t', encoding='utf-8')

    # else:
    #     # remove config file XML
    #     # if (RemoveConfigFile): os.remove( pathlib.Path(ConfigFile) )
    #     # write the stats in a file and remove the folder
    #     if (SummaryFunction):
    #         OutputFolder = self.get_outputPath(SampleID, ReplicateID)
    #         SummaryFile = self.outputs_folder+'SummaryFile_%06d_%02d.csv'%(SampleID,ReplicateID)
    #         ParamNames = [self.parameters[param_key][1] for param_key in self.keys_variable_params]
    #         dic_params = {ParamNames[i]: Parameters[i] for i in range(len(Parameters))}
    #         try: result_summary = SummaryFunction(OutputFolder,SummaryFile, dic_params,  SampleID, ReplicateID)
    #         except OSError as error: print(f"\t{error}\n\tError in SummaryFunction! (Sample: {SampleID} and Replicate: {ReplicateID}).")
    #         except SystemError as error: print(f"\t{error}\n\tError in SummaryFunction! (Sample: {SampleID} and Replicate: {ReplicateID}).")
    # return 0

if __name__ == '__main__':
    # PhysiCell_Model, Samples, Replicates = args_run_simulations(sys.argv[1:]) 

    Samples = generate_list_of_config_files(sys.argv[1])
    # make_output_directories(Samples) --> in create_files.py now

    NumSimulations = len(Samples) # - which is a list of strings number of samples - this is in fact the number of simulations your run - 
            # ID of each aparmaeter set. Would need the right folder for storage and EVERYTHING is ready to run. And the physicell executable is in the right place.
            # also need to chagne the seed. 
            # I THINK samples is indexless. 
    NumSimulationsPerRank  = int(NumSimulations/size)

    # print(NumSimulationsPerRank)
    
    print(f"Total number of simulations: {NumSimulations} Simulations per rank: {NumSimulationsPerRank}")

    data = np.array_split(np.arange(NumSimulations),size, axis=0) # split [0,1,...,NumSimulations-1] in size arrays equally (or +1 in some ranks)

    print(data)

    ## We won't hae "replicates" in the sense that its here - this will all be managed prior to running by number of config files generated

    for ind_sim in data[rank]:
        sampleID = Samples[ind_sim]
        print(sampleID)
        # replicateID = Replicates[ind_sim]
        # print('Rank: ',rank, ', Simulation: ', ind_sim, ', Sample: ', sampleID,', Replicate: ', replicateID)
        print('Rank: ',rank, ', Simulation: ', ind_sim, ', Sample: ', sampleID)
        # model(PhysiCell_Model.get_configFilePath(sampleID, replicateID), PhysiCell_Model.executable)
        # in the original setup, this (ths python file) would have been called with several arguments. Here, we only need to call it with the directory with the config
        # files in them. And the executable is currently hard coded. 
        model(sampleID, './PhysiBoSS_Cell_Lines')
