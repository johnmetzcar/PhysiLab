# This is all (the tools, the base notebook, the techniques) from the PhysiBoSS Tutorial paper:

Ruscone, M., Checcoli, A., Heiland, R., Barillot, E., Macklin, P., Calzone, L., Noël, V., n.d. Building multi-scale models with PhysiBoSS, an agent-based modeling tool.

https://github.com/PhysiBoSS/PhysiBoSS/blob/development/sample_projects_intracellular/boolean/tutorial/paper/PhysiBoSS_tutorial_main_text.pdf

Modified for the PhysiLab project by John Metzcar

---------------------------------- Here are the original instructions ---------------------------------------------------

# PhysiBoSS Tutorial Jupyter Notebook 

In this folder we provide a few jupyter notebooks which analyse either the intracellular Boolean models or the PhysiBoSS models. 

## Available notebooks

CellFate_Analysis.ipynb : Analysis of the Calzone's cell fate model. 

Cell_cycle_boolean_analysis.ipynb : Analysis of the results of the cell cycle intracellular model

Cell_cycle_analysis.ipynb : Analysis of the results of the cell cycle PhysiBoSS model


## Running the notebooks

To install dependencies, you can use conda and run :

`conda env create --file environment.yml`

and then activate your environment with :

`conda activate tutorial-notebooks`

Finally, you can run Jupyter Lab and open the notebooks :

`jupyter-lab`


