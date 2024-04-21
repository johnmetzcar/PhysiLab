# for help with ElementTree: https://docs.python.org/3/library/xml.etree.elementtree.html
# creates maboss (.bnd) AND PhysiCell (.xml) files for controls imported from csv files
# Work led by Katie Pletz with contributions from John Metzcar
# Katie Pletz and John Metzcar 2023-2024. 


import os
import os.path
import sys
import pandas as pd
import xml.etree.ElementTree as ET

###################### EDIT THESE VARIABLES ######################
### Note - eventually will put this in a config file (.ini) or something similar

######################## Base Model Files ########################

# MaBoSS files
# EDIT TO BE NAME OF ORIGINAL BND FILE (MaBoSS model)
MaBoSS_file = "TLGL_base.bnd"
# EDIT TO BE NAME OF ORIGINAL CFG FILE (MaBoSS model)
MaBoSS_config = "TLGL_base_survival_attractors.cfg"

# PhysiCell files
# EDIT TO BE NAME OF ORIGINAL XML FILE (PhysiCell model)
PhysiCell_file = "base_spatial_model_file_multiple_interventions.xml"
print("CURRENTLY OTHER PHYSICELL CONFIG FILES ARE ASSUMED TO BE IN CONFIG (rules, cell initialization, Dirichelet nodes, and substrate initializaiton)" )

#################### Intervention Files ########################

stable_motifs_interventions_file = 'FormattedInternalMergeResults.csv'
ibmfa_interventions_file = 'IBMFA_top_interventions.csv'
edgetic_interventions_file = 'single_edge_perturbations_top_interventions.csv'

################### Pathing Variables ########################

# relative path to input directory
# EDIT TO BE LOCATION OF ORIGINAL MODEL FILES
rel_input_dir = 'Variant_Model_Files/'

# relative path to output directory
# EDIT TO BE LOCATION OF WHERE YOU WANT TO STORE THE NEW MODEL FILES
rel_output_dir = 'leukemia_spatial_model_files/'

rel_simulation_output_dir = 'leukemia_spatial_output/'

########################### PATHING ###########################
# move to PhysiCell root directory - this assumes this script is one folder below the root
os.chdir("../")
full_path = os.getcwd()
print(full_path)

# INPUTS
input_dir = os.path.join(full_path, rel_input_dir)
print("Input directory: ")
print(input_dir)
input_MaBoSS_file = os.path.join(input_dir, MaBoSS_file)
input_MaBoSS_config = os.path.join(input_dir, MaBoSS_config)

# MODELS
print('Base model files:')
print("Input BND file: ")
print(input_MaBoSS_file)
print("Input CFG file: ")
print(input_MaBoSS_config)

input_PC_file = os.path.join(input_dir, PhysiCell_file)
print("Input XML file: ")
print(input_PC_file)

# Interventions 
print("Interventions Files:")
stable_motifs_file_and_path = os.path.join(input_dir, stable_motifs_interventions_file)
print("Stable Motifs File: ")
print(stable_motifs_file_and_path)

print("IBMFA File: ")
ibmfa_file_and_path = os.path.join(input_dir, ibmfa_interventions_file)
print(ibmfa_file_and_path)

print("Edgetic File: ")
edgetic_file_and_path = os.path.join(input_dir, edgetic_interventions_file)
print(edgetic_file_and_path)

# OUTPUTS
output_dir = os.path.join(full_path, rel_output_dir)
print("Output directory: ")
print(output_dir)

simulation_output_dir = os.path.join(full_path, rel_simulation_output_dir)
print("Simulation output directory: ")
print(simulation_output_dir)

input("Press Enter to continue...\n Press Ctrl + C to exit...")

# make sure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# make sure the simulation output directory exists
if not os.path.exists(simulation_output_dir):
    os.makedirs(simulation_output_dir)
    
# move to the output directory
os.chdir(output_dir)


# function for reading bnd file
def getBNDdata():
    # open original bnd  --> location and name set in INPUT section above
    baseBND = open(input_MaBoSS_file, "r+")
    bndText = baseBND.readlines()

    # create one big long string of the entire bnd file
    bndString = ""
    for j in range(len(bndText)):
        bndString = bndString + bndText[j]

    # separate the string into a list of strings: one for each node
    nodes = bndString.split("Node")

    # remove the empty string at the beginning
    nodes.pop(0)

    # create a dictionary with key = node and value = logic
    nodeDict = {}
    for node in nodes:
        nodename = node[1:].split(" {")[0]
        nodelogic = node.split("\n")[1]
        nodelogic = nodelogic.split(";")[0]
        nodeDict[nodename] = nodelogic

    return nodeDict

def getCFG():
    # open original cfg  --> location and name set in INPUT section above
    baseCFG = open(input_MaBoSS_config, "r+")
    cfgText = baseCFG.readlines()

    return cfgText


## May eventually have this function take xml file as input and modify it
def createXML(intervention, substrateNames, numInterventions, decayID, replicateID):
    # THIS FUNCTION ALSO MAKES THE OUTPUT DIRECTORIES!!!!!!!!!
    # lets add doc strings eventually ... 

    # load and parse base xml file --> location and name set in INPUT section above
    base_xml = open(input_PC_file)
    tree = ET.parse(base_xml)
    xml_root = tree.getroot()

    # modify save
    # use first intervention
    save = xml_root.find("save")
    SVG = save.find("SVG")
    plot_substrate = SVG.find("plot_substrate")
    substrate = plot_substrate.find("substrate")
    substrate.text = substrateNames[0]      

    # microenvironment variable names
    microenvironment = xml_root.find("microenvironment_setup")
    variable1 = microenvironment.findall("variable")[0]
    variable2 = microenvironment.findall("variable")[1]
    variable3 = microenvironment.findall("variable")[2]

    variable1.attrib = {"name": substrateNames[0], "units": "dimensonless", "ID": "0"}
    physical_parameter_set = variable1.find("physical_parameter_set")
    decay_rate = physical_parameter_set.find("decay_rate")
    if decayID == "1":
        decay_rate.text = "0.00475" 
    # elif decayID == "2":
    #     decay_rate.text = "0.00192" # 6 hours above 0.5 - current threshold 
    # else:
    #     decay_rate.text = "0.000963" # 12 hours above 0.5 - current threshold

    if numInterventions > 1:
        variable2.attrib = {"name": substrateNames[1], "units": "dimensionless", "ID": "1"}
        physical_parameter_set = variable2.find("physical_parameter_set")
        decay_rate = physical_parameter_set.find("decay_rate")
        if decayID == "1":
            decay_rate.text = "0.00475" 
            # elif decayID == "2":
            #     decay_rate.text = "0.00192" # 6 hours above 0.5 - current threshold 
            # else:
            #     decay_rate.text = "0.01" "0.000963" # 12 hours above 0.5 - current threshold
    else:
        microenvironment.remove(variable2)

    if numInterventions == 3:
        variable3.attrib = {"name": substrateNames[2], "units": "dimensionless", "ID": "2"}
        physical_parameter_set = variable3.find("physical_parameter_set")
        decay_rate = physical_parameter_set.find("decay_rate")
        if decayID == "1":
            decay_rate.text = "0.00475" 
        # elif decayID == "2":
        #     decay_rate.text = "0.00192" # 6 hours above 0.5 - current threshold 
        # else:
        #     decay_rate.text = "0.01" "0.000963" # 12 hours above 0.5 - current threshold
    else:
        microenvironment.remove(variable3)

    # modify chemotaxis
    cell_defs = xml_root.find("cell_definitions")
    cell_def = cell_defs.find("cell_definition")
    phenotype = cell_def.find("phenotype")
    motility = phenotype.find("motility")
    options = motility.find("options")
    chemotaxis = options.find("chemotaxis")
    substrate = chemotaxis.find("substrate")
    substrate.text = substrateNames[0]

    # modify advanced chemotaxis
    adv_chemotaxis = options.find("advanced_chemotaxis")
    sensitivities = adv_chemotaxis.find("chemotactic_sensitivities")
    sensitivity = sensitivities.findall("chemotactic_sensitivity")
    sensitivity1 = sensitivity[0]
    sensitivity2 = sensitivity[1]
    sensitivity3 = sensitivity[2]

    sensitivity1.attrib = {"substrate": substrateNames[0]}
    if numInterventions > 1:
        sensitivity2.attrib = {"substrate": substrateNames[1]}
    else:
        sensitivities.remove(sensitivity2)
    
    if numInterventions == 3:
        sensitivity3.attrib = {"substrate": substrateNames[2]}
    else:
        sensitivities.remove(sensitivity3)

    # modify secretion
    secretion = phenotype.find("secretion")
    substrates = secretion.findall("substrate")
    substrate1 = substrates[0]
    substrate2 = substrates[1]
    substrate3 = substrates[2]

    substrate1.attrib = {"name": substrateNames[0]}
    if numInterventions > 1:
        substrate2.attrib = {"name": substrateNames[1]}
    else:
        secretion.remove(substrate2)
    
    if numInterventions == 3:
        substrate3.attrib = {"name": substrateNames[2]}
    else:
        secretion.remove(substrate3)

    # change filenames
    intracellular = phenotype.find("intracellular")
    bndfile = intracellular.find("bnd_filename")
    
    bndfile.text = rel_output_dir + intervention + ".bnd"
    cfgfile = intracellular.find("cfg_filename")
    cfgfile.text = rel_output_dir  + intervention + ".cfg"

    # modify mapping info
    mapping = intracellular.find("mapping")
    inputs = mapping.findall("input")
    input1 = inputs[0]
    input2 = inputs[1]
    input3 = inputs[2]
    
    input1.attrib = {"physicell_name": substrateNames[0], "intracellular_name": substrateNames[0]}
    settings = input1.find("settings")
    action = settings.find("action")
    action.text = "activation"

    if numInterventions > 1:
        input2.attrib = {"physicell_name": substrateNames[1], "intracellular_name": substrateNames[1]}
        settings = input2.find("settings")
        action = settings.find("action")
        action.text = "activation"
    else:
        # remove input 2
        mapping.remove(input2)
    
    if numInterventions == 3:
        input3.attrib = {"physicell_name": substrateNames[2], "intracellular_name": substrateNames[2]}
        settings = input3.find("settings")
        action = settings.find("action")
        action.text = "activation"
    else:
        # remove substrate 3
        mapping.remove(input3)
    
    # change user parameters
    user_params = xml_root.find("user_parameters")
    substrate1 = user_params.find("substrate_name1")
    substrate2 = user_params.find("substrate_name2")
    substrate3 = user_params.find("substrate_name3")
    substrate1.text = substrateNames[0]

    if numInterventions > 1:
        substrate2.text = substrateNames[1]
    else:
        # remove substrate 2
        user_params.remove(substrate2)
    
    if numInterventions == 3:
        substrate3.text = substrateNames[2]
    else:
        # remove substrate 3
        user_params.remove(substrate3)
    
    num_substrates = user_params.find("num_substrates")
    num_substrates.text = str(numInterventions)

    # create serial number
    modelID = decayID + replicateID

    # modify output file
    folder = save.find("folder")
    folder.text = rel_simulation_output_dir + intervention + "_" + modelID

    # change seed value
    random_seed = user_params.find("random_seed")
    random_seed.text = modelID

    # create new xml file
    # output directory given above in the OUTPUT section (unless changed in MAIN script)
    fileName = intervention + "_" + modelID + ".xml"
    tree.write(fileName)
    print("Created file " + fileName)

    # make the output directory
    output_dir = os.path.join(simulation_output_dir, intervention + "_" + modelID)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def createCFG(intervention, substrateNames):
    # part1 = ["PDGF.istate = 0;\n", 
    #             "IL15.istate = 1;\n", 
    #             "Stimuli.istate = 1;\n",
    #             "Stimuli2.istate = 0;\n",
    #             "CD45.istate = 0;\n",
    #             "TAX.istate = 0;\n"]

    # part2 = ["Apoptosis.istate = 0;\n",
    #             "discrete_time = 0;\n",
    #             "use_physrandgen = FALSE;\n",
    #             "// seed_pseudorandom = 100;\n",
    #             "sample_count = 100000;\n\n",
    #             "max_time = 2000.0;\n",
    #             "time_tick = 1.0;\n\n",
    #             "thread_count = 8;"]

    part1 = getCFG()

    cfgName = intervention + ".cfg"
    # output directory given above in the OUTPUT section (unless changed in MAIN script)
    cfgfile = open(cfgName, "w")
    cfgfile.writelines(part1)
    for s in substrateNames:
        newNode = s + ".istate = 0;\n"
        cfgfile.write(newNode)
    # cfgfile.writelines(part2)

    print("Created file " + intervention + ".cfg")


###################################################
############## END FUNCTIONS ######################
    
############# START OF MAIN SCRIPT ################
    
# import interventions files
ibmfa = pd.read_csv(ibmfa_file_and_path)
stableMotifs = pd.read_csv(stable_motifs_file_and_path, header = None)
edgetic = pd.read_csv(edgetic_file_and_path)

###################################################
############## IBMFA Interventions ################
###################################################

# if you want to store the IBMFA files in the same directory as the other files, comment out the next lines (NOT TESTED)
    
# rel_output_dir_IBMFA = 'leukemia_model_files/IBMFA_files/'
# output_dir_IBMFA = os.path.join(full_path, rel_output_dir_IBMFA)
# print("Output directory: ")
# print(output_dir_IBMFA)
# os.chdir(output_dir_IBMFA)

ibmfaDict = {}
for i in range(len(ibmfa)):
    # create dictionary entries with keys = intervention, values = subinterventions
    intervention = (ibmfa.iloc[i, 0]).replace(" ", "")
    subinterventions = intervention.split("&")
    ibmfaDict[intervention] = subinterventions

for intervention in ibmfaDict:
    nodeDict = getBNDdata()
    numInterventions = len(ibmfaDict[intervention])
    substrateNames = []
    interventionName = ""
    for subintervention in ibmfaDict[intervention]:
        nodeOfInterest = subintervention.split("-")[0]
        statusOfInterest = subintervention.split("-")[1]

        if(statusOfInterest == "0"):
            substrate = "anti_" + nodeOfInterest
        else:
            substrate = "pro_" + nodeOfInterest
        substrateNames.append(substrate)
        interventionName += substrate

        newLogic = nodeDict.get(nodeOfInterest) + ")"
        if (statusOfInterest == "0"):
            # suppression
            physiName = "anti_" + nodeOfInterest
            newLogic = newLogic[:-1] + " & !" + physiName
        else:
            # promotion
            physiName = "pro_" + nodeOfInterest
            newLogic = newLogic[:-1] + " | " + physiName
        
        # edit logic for node of interest
        nodeDict[physiName] = "  logic = " + physiName
        nodeDict[nodeOfInterest] = newLogic

    ibmfaDict[intervention] = substrateNames

    # add algorithm name to intervention name
    interventionName = "IB_" + interventionName

    # create bnd file
    fileName = interventionName + ".bnd"
    newFile = open(fileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? 1.0 : 0;\n  rate_down = @logic ? 0 : 1.0 ;\n}"
        newFile.write(nodeString + "\n\n")
    print("Created file " + fileName)

    # create cfg file
    createCFG(interventionName, ibmfaDict[intervention])

    # create xml files
    # for d in range(3): # makes 3 decay rates
    #     for r in range(3): # makes 3 replicates
    #         decay = str(d + 1)
    #         replicate = str(r + 1)
            # createXML(interventionName, substrateNames, numInterventions, decay, replicate)
    createXML(interventionName, substrateNames, numInterventions, "1", "1")

print("--------------------------")
print("Finished with IBMFA files.")
print("--------------------------")

###################################################
############ Stable Motif Interventions ############
###################################################



# if you want to store the SM files in the same directory as the other files, comment out the next lines (NOT TESTED)
    
# rel_output_dir_SM = 'leukemia_model_files/SM_files/'
# output_dir_SM = os.path.join(full_path, rel_output_dir_SM)
# print("Output directory: ")
# print(output_dir_SM)
# os.chdir(output_dir_SM)

## Can this go into a function???
# split into interventions and subinterventions
stableMotifDict = {}
for i in range(len(stableMotifs)):
    # create dictionary entries with keys = intervention, values = subinterventions
    intervention = (stableMotifs.iloc[i, 0]).replace(" ", "")
    subinterventions = intervention.split("&")
    stableMotifDict[intervention] = subinterventions

for intervention in stableMotifDict:
    nodeDict = getBNDdata()
    numInterventions = len(stableMotifDict[intervention])
    substrateNames = []
    interventionName = ""
    for subintervention in stableMotifDict[intervention]:
        nodeOfInterest = subintervention.split("-")[0]
        statusOfInterest = subintervention.split("-")[1]

        if(statusOfInterest == "0"):
            substrate = "anti_" + nodeOfInterest
        else:
            substrate = "pro_" + nodeOfInterest
        substrateNames.append(substrate)
        interventionName += substrate
    
        # edit bnd file
        newLogic = nodeDict.get(nodeOfInterest) + ")"
        if (statusOfInterest == "0"):
            # suppression
            physiName = "anti_" + nodeOfInterest
            newLogic = newLogic[:-1] + " & !" + physiName
        else:
            # promotion
            physiName = "pro_" + nodeOfInterest
            newLogic = newLogic[:-1] + " | " + physiName

        # update logic for node of interest
        nodeDict[physiName] = "  logic = " + physiName
        nodeDict[nodeOfInterest] = newLogic

    stableMotifDict[intervention] = substrateNames
    
    # add algorithm name to intervention name
    interventionName = "SM_" + interventionName

    # create bnd file
    # for each subintervention, add the corresponding input node and update logic
    fileName = interventionName + ".bnd"
    newFile = open(fileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? 1.0 : 0;\n  rate_down = @logic ? 0 : 1.0 ;\n}"
        newFile.write(nodeString + "\n\n")
    print("Created file " + fileName)

    # create cfg file
    createCFG(interventionName, stableMotifDict[intervention])

    # create xml files
    # for d in range(3): # makes 3 decay rates
    #     for r in range(3): # makes 3 replicates
    #         decay = str(d + 1)
    #         replicate = str(r + 1)
            # createXML(interventionName, substrateNames, numInterventions, decay, replicate)
    createXML(interventionName, substrateNames, numInterventions, "1", "1")


print("----------------------------------")
print("Finished with Stable Motifs files.")
print("----------------------------------")


###################################################
############ Single Edge Perturbations ############
###################################################

# create bnd files for edgetic perturbations results

# if you want to store the edgetic perturbation files in the same directory as the other files, comment out the next lines (NOT TESTED)
    
# rel_output_dir_EG = 'leukemia_model_files/EG_files/'
# output_dir_EG = os.path.join(full_path, rel_output_dir_EG)
# print("Output directory: ")
# print(output_dir)
# os.chdir(output_dir_EG)

for i in range(len(edgetic)):
    # separate components of result
    source = edgetic.iloc[i, 0].split("-->")[0]
    target = (edgetic.iloc[i, 0].split("-->")[1]).split(",")[0]
    action = (edgetic.iloc[i, 0].split("-->")[1]).split(",")[3]

    # get base bnd file info
    nodeDict = getBNDdata()

    # update logic of nodes for intervention
    # update existing logic for target node
    if (action == "0"):
        # suppression
        newName = "anti_" + source + "_" + target
        nodeDict[newName] = "  logic = " + newName
        #nodeDict[target] = nodeDict[target] + " & !" + newName
        oldLogic = nodeDict[target]
        modified = "(" + source + " & !" + newName + ")"
        newLogic = oldLogic.replace(source, modified)
        nodeDict[target] = newLogic
    else:
        # excitation
        newName = "pro_" + source + "_" + target
        nodeDict[newName] = "  logic = " + newName
        #nodeDict[target] = nodeDict[target] + " | " + newName
        oldLogic = nodeDict[target]
        modified = "(" + source + " | " + newName + ")"
        newLogic = oldLogic.replace(source, modified)
        nodeDict[target] = newLogic
    
    fileName = source + "_" + target + "_" + action

    # add algorithm name to intervention name
    fileName = "EG_" + fileName # differs slightly from previous names - may need to change

    # create new xml files
    # def createXML(intervention, substrateNames, numInterventions, decayID, replicateID):
    # for d in range(3):
    #     for r in range(3):
    #         decay = str(d + 1)
    #         replicate = str(r + 1)
    #         createXML(fileName, [newName], 1, decay, replicate)
    createXML(fileName, [newName], 1, "1", "1")
            
    # create new BND file
    newFileName = fileName + ".bnd"
    newFile = open(newFileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? 1.0 : 0;\n  rate_down = @logic ? 0 : 1.0 ;\n}"
        newFile.write(nodeString + "\n\n")

    print("Created file " + newFileName)

    # create cfg file
    createCFG(fileName, [newName])
    print("Created file " + fileName + ".cfg")


print("----------------------------------------------")
print("Finished with Single Edge Perturbations Files.")
print("----------------------------------------------")

