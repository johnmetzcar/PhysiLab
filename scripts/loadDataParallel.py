## Loads population data from PhysiCell outputs and exports to csv
## Runs quickly did 250 directories in less than 20 minutes on Big Red login node. 
## Would be nice to have it take CLI arguments - input and output directories and final file name.
# BE CAREFUL WITH USE - could overwrite files. 
## Copy write Katie Pletz and John Metzcar 2024.


import os
import multiprocessing
import pcdl
import pandas as pd
import sys

def loadData(args):
    intervention, input_dir, output_dir = args
    intervention_input_path = os.path.join(input_dir, intervention)

    print("Processing:", intervention)
    print("Input Path:", intervention_input_path)

    mcds_ts = pcdl.TimeSeries(intervention_input_path, microenv=False, settingxml=None, graph=False, verbose=False)

    data = {'intervention': intervention}
    for mcds in mcds_ts.get_mcds_list():
        df_cell = mcds.get_cell_df()
        live_cells = len(df_cell[df_cell['dead'] == False])
        data[mcds.get_time()] = live_cells

    df = pd.DataFrame([data])
    output_csv_path = os.path.join(intervention_input_path, 'live_cells.csv')
    df.to_csv(output_csv_path, index=False)

def aggregateResults(input_dir, output_dir, interventions, output_file_name):
    aggregated_df = pd.DataFrame()

    for intervention in interventions:
        csv_path = os.path.join(input_dir, intervention, 'live_cells.csv')
        if os.path.exists(csv_path):
            data = pd.read_csv(csv_path)
            aggregated_df = pd.concat([aggregated_df, data], ignore_index=True)

    output_csv_path = os.path.join(output_dir, output_file_name)
    aggregated_df.to_csv(output_csv_path, index=False)
    print("Aggregated results saved to:", {output_csv_path})

def main(input_dir, output_dir, output_file_name):
    # Assuming the below directory adjustments and path listings are similar to your initial setup
    rel_input_dir = input_dir
    rel_output_dir = output_dir

    full_path = os.getcwd()
    print("Current Working Directory:", full_path)

    input_dir = os.path.join(full_path, rel_input_dir)
    print("Input Directory:", input_dir)

    output_dir = os.path.join(full_path, rel_output_dir)
    print("Output Directory:", output_dir)

    # input("Press Enter to continue...\n Press Ctrl + C to exit...")

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    interventions = [f.name for f in os.scandir(input_dir) if f.is_dir()]
    interventions = sorted(interventions)

    args_list = [(intervention, input_dir, output_dir) for intervention in interventions]

    # Use multiprocessing to process each intervention in parallel
    with multiprocessing.Pool() as pool:
        pool.map(loadData, args_list)

    # Aggregate results after all processes have completed
    aggregateResults(input_dir, output_dir, interventions, output_file_name)

if __name__ == '__main__':
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    output_file_name = sys.argv[3]
    main(input_dir, output_dir, output_file_name)