# Creates grid of population plots organized by intervention for a user-defined subset of interventions

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# location of population dataframe
df = pd.read_csv("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\PhysiLab Test Copy\\Dataframes\\live_cells.csv")

# get list of all unique interventions
interventions = df['intervention'].values.tolist()

# get list of unique "parent" interventions
df['parent_intervention'] = df['intervention'].str[0:-3]
parents = df.parent_intervention.unique().tolist()

# all possible serial numbers
kids = ["11", "12", "13", "21", "22", "23", "31", "32", "33"]

# define dimensions for grid (rows, cols)
figure, axis = plt.subplots(9, 7)

# change range to plot desired interventions:
# lower bound = index of first intervention (starting from 0)
# upper bound = index of last intervention + 1
# e.g., range(7, 14) plots interventions numbered 8-14 (indices 7-13)
for p in range(7, 14):
    # get all interventions with matching parent intervention
    children = df.loc[df['parent_intervention'] == parents[p]]
    children = (children['intervention']).values.tolist()

    # for each child
    for i in range(len(children)):   
        # get population data and convert to numpy arrays
        data = df.loc[df['intervention'] == children[i]]
        data = data.loc[:, data.columns != 'intervention']
        data = data.loc[:, data.columns != 'parent_intervention']
        data = data.loc[:, data.columns != 'Unnamed: 0']
        time = np.array(data.columns)
        population = np.array(data.values.flatten().tolist())

        # get row and column indices
        row = i
        col = p - 7

        # create plot
        axis[row, col].plot(time, population)
        if row == 0:
            axis[row, col].set_title(parents[p])
        if col == 0:
            axis[row, col].set_ylabel(kids[i])
        axis[row, col].set_xticklabels([])
        axis[row, col].set_xticks([])
        axis[row, col].set_yticklabels([])
        axis[row, col].set_yticks([])
    
    print("Adding plot " + children[i])

plt.show()