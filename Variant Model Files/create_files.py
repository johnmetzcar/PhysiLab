# creates maboss (.bnd) AND PhysiCell (.xml) files for controls imported from csv files

import os
import sys
import pandas
import xml.etree.ElementTree as ET

# function for reading bnd file
def getBNDdata():
    # open original bnd  --> EDIT TO BE LOCATION OF ORIGINAL BND FILE
    baseBND = open("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\TLGL.bnd", "r+")
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

def createXML(intervention, substrateNames, numInterventions, decayID, replicateID):
    # load and parse base xml file
    base_xml = open("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\base_model_file_two_interventions.xml")
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
        decay_rate.text = "1"
    elif decayID == "2":
        decay_rate.text = "0.1"
    else:
        decay_rate.text = "0.01"

    if numInterventions > 1:
        variable2.attrib = {"name": substrateNames[1], "units": "dimensionless", "ID": "0"}
        physical_parameter_set = variable2.find("physical_parameter_set")
        decay_rate = physical_parameter_set.find("decay_rate")
        if decayID == "1":
            decay_rate.text = "1"
        elif decayID == "2":
            decay_rate.text = "0.1"
        else:
            decay_rate.text = "0.01"
    else:
        microenvironment.remove(variable2)

    if numInterventions == 3:
        variable3.attrib = {"name": substrateNames[2], "units": "dimensionless", "ID": "0"}
        physical_parameter_set = variable3.find("physical_parameter_set")
        decay_rate = physical_parameter_set.find("decay_rate")
        if decayID == "1":
            decay_rate.text = "1"
        elif decayID == "2":
            decay_rate.text = "0.1"
        else:
            decay_rate.text = "0.01"
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
    bndfile.text = "config/" + intervention + ".bnd"
    cfgfile = intracellular.find("cfg_filename")
    cfgfile.text = "config/" + intervention + ".cfg"

    # modify mapping info
    mapping = intracellular.find("mapping")
    inputs = mapping.findall("input")
    input1 = inputs[0]
    input2 = inputs[1]
    input3 = inputs[2]
    
    input1.attrib = {"physicell_name": substrateNames[0], "intracellular_name": substrateNames[0]}
    settings = input1.find("settings")
    action = settings.find("action")
    if substrateNames[0].startswith("pro"):
        action.text = "activation"
    else:
        action.text = "inhibition"

    if numInterventions > 1:
        input2.attrib = {"physicell_name": substrateNames[1], "intracellular_name": substrateNames[1]}
        settings = input2.find("settings")
        action = settings.find("action")
        if substrateNames[1].startswith("pro"):
            action.text = "activation"
        else:
            action.text = "inhibition"
    else:
        # remove input 2
        mapping.remove(input2)
    
    if numInterventions == 3:
        input3.attrib = {"physicell_name": substrateNames[2], "intracellular_name": substrateNames[2]}
        settings = input3.find("settings")
        action = settings.find("action")
        if substrateNames[2].startswith("pro"):
            action.text = "activation"
        else:
            action.text = "inhibition"
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
    folder.text = intervention + "_" + modelID

    # change seed value
    random_seed = user_params.find("random_seed")
    random_seed.text = modelID

    # create new xml file
    fileName = intervention + "_" + modelID + ".xml"
    tree.write(fileName)
    print("Created file " + fileName)

def createCFG(intervention, substrateNames):
    part1 = ["PDGF.istate = 0;\n", 
                "IL15.istate = 1;\n", 
                "Stimuli.istate = 1;\n",
                "Stimuli2.istate = 0;\n",
                "CD45.istate = 0;\n",
                "TAX.istate = 0;\n"]

    part2 = ["Apoptosis.istate = 0;\n",
                "Proliferation.istate = 0;\n\n",
                "$time_scale = 0.1;\n",
                "$oxygen_concentration = 0.0;\n",
                "$glucose_concentration = 0.0;\n\n",
                "discrete_time = 1;\n",
                "use_physrandgen = FALSE;\n",
                "// seed_pseudorandom = 100;\n",
                "sample_count = 1;\n\n",
                "max_time = 1.0;\n",
                "time_tick = 0.0004;\n\n",
                "thread_count = 1;"]

    cfgName = intervention + ".cfg"
    cfgfile = open(cfgName, "w")
    cfgfile.writelines(part1)
    for s in substrateNames:
        newNode = s + ".istate = 0;\n"
        cfgfile.write(newNode)
    cfgfile.writelines(part2)

    print("Created file " + intervention + ".cfg")

# set working directory --> EDIT TO BE LOCATION OF YOUR CSV FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files")

# import csv files
ibmfa = pandas.read_csv("IBMFA_top_interventions.csv")
stableMotifs = pandas.read_csv("FormattedInternalMergeResults.csv", header = None)
edgetic = pandas.read_csv("single_edge_perturbations_top_interventions.csv")

# create files for IBMFA results

# change directory --> EDIT TO BE THE LOCATION YOU WANT TO STORE IBMFA FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\IBMFA Files")

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

    # create bnd file
    fileName = interventionName + ".bnd"
    newFile = open(fileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")
    print("Created file " + fileName)

    # create cfg file
    createCFG(interventionName, ibmfaDict[intervention])

    # create xml files
    for d in range(3):
        for r in range(3):
            decay = str(d + 1)
            replicate = str(r + 1)
            createXML(interventionName, substrateNames, numInterventions, decay, replicate)

print("--------------------------")
print("Finished with IBMFA files.")
print("--------------------------")

# create files for Stable Motifs results

# change directory --> EDIT TO BE THE LOCATION YOU WANT TO STORE STABLE MOTIFS FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\Stable Motifs Files")

stableMotifDict = {}
for i in range(len(stableMotifs)):
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

    # create bnd file
    fileName = interventionName + ".bnd"
    newFile = open(fileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")
    print("Created file " + fileName)

    # create cfg file
    createCFG(interventionName, stableMotifDict[intervention])

    # create xml files
    for d in range(3):
        for r in range(3):
            decay = str(d + 1)
            replicate = str(r + 1)
            createXML(interventionName, substrateNames, numInterventions, decay, replicate)

print("----------------------------------")
print("Finished with Stable Motifs files.")
print("----------------------------------")

# create bnd files for edgetic perturbations results

# change directory --> EDIT TO BE THE PLACE YOU WANT TO STORE SINGLE EDGE PERTURBATION FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\Single Edge Perturbations Files")

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

    # create new xml files
    # def createXML(intervention, substrateNames, numInterventions, decayID, replicateID):
    for d in range(3):
        for r in range(3):
            decay = str(d + 1)
            replicate = str(r + 1)
            #createXML(fileName, statusOfInterest, decay, replicate)
            createXML(fileName, [newName], 1, decay, replicate)

    # create new BND file
    newFileName = fileName + ".bnd"
    newFile = open(newFileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")

    print("Created file " + newFileName)

    # create cfg file
    createCFG(fileName, [newName])
    print("Created file " + fileName + ".cfg")

print("----------------------------------------------")
print("Finished with Single Edge Perturbations Files.")
print("----------------------------------------------")
