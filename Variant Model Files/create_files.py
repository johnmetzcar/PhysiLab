# creates maboss (.bnd) AND PhysiCell (.xml and .cfg) files for controls imported from csv files

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

# function for adding new layers for xml file
# doesn't return anything, just modifies the tree data
def addMappingInfo(root, status, physicellName):

    # locate intracellular layer
    cell_defs = root.find("cell_definitions")
    cell_def = cell_defs.find("cell_definition")
    phenotype = cell_def.find("phenotype")
    intracellular = phenotype.find("intracellular")
    mapping = intracellular.find("mapping")
    input = mapping.find("input")

    # modify physicell_name and intracellular_name
    input.attrib = {"physicell_name": physicellName, "intracellular_name": physicellName}

    # modify action
    settings = input.find("settings")
    action = settings.find("action")
    if status == "0":
        action.text = "inhibition"
    else:
        action.text = "activation"

    # add mapping info
    #mapping = ET.SubElement(intracellular, "mapping", {})
    #input = ET.SubElement(mapping, "input", {"physicell_name": physicellName, "intracellular_name": physicellName})
    #settings = ET.SubElement(input, "settings", {})
    #action = ET.SubElement(settings, "action", {})
    #if status == "0":
    #    action.text = "inhibition"
    #else:
    #    action.text = "activation"
    #threshold = ET.SubElement(settings, "threshold", {})
    #threshold.text = "0.5"
    #smoothing = ET.SubElement(settings, "smoothing", {})
    #smoothing.text = "0"

def createCFG(fileName):
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

    # new node setting
    newNode = fileName + ".istate = 0;\n"

    cfgName = fileName + ".cfg"
    cfgfile = open(cfgName, "w")
    cfgfile.writelines(part1)
    cfgfile.write(newNode)
    cfgfile.writelines(part2)

# set working directory --> EDIT TO BE LOCATION OF YOUR CSV FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files")

# import csv files
ibmfa = pandas.read_csv("IBMFA_top_interventions.csv")
stableMotifs = pandas.read_csv("FormattedInternalMergeResults.csv", header = None)
edgetic = pandas.read_csv("single_edge_perturbations_top_interventions.csv")

# create files for IBMFA results

# change directory --> EDIT TO BE THE LOCATION YOU WANT TO STORE IBMFA FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\IBMFA Files")

for i in range(len(ibmfa)):
    # get base bnd file info
    nodeDict = getBNDdata()

    # load and parse base xml file
    base_xml = open("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\base_model_file.xml")
    tree = ET.parse(base_xml)
    xml_root = tree.getroot()

    # get node(s) of interest and desired state (from csv)
    intervention = (ibmfa.iloc[i, 0]).replace(" ", "")
    subinterventions = intervention.split("&")
    fileName = ""

    for s in subinterventions:
        nodeOfInterest = s.split("-")[0]
        statusOfInterest = s.split("-")[1]

        # create new node for intervention
        newName = "pro_" + nodeOfInterest
        nodeDict[newName] = "  logic = " + newName  

        # find node logic corresponding to node of interest
        newLogic = nodeDict.get(nodeOfInterest) + ")"

        if(statusOfInterest == "0"):
            # suppression
            newLogic = newLogic[:-1] + " & !" + newName
            physiName = "anti_" + nodeOfInterest
            fileName += physiName
        else:
            # promotion
            newLogic = newLogic[:-1] + " | " + newName
            physiName = "pro_" + nodeOfInterest
            fileName += physiName

        # edit logic for node of interest
        nodeDict[nodeOfInterest] = newLogic

    # add xml mapping data
    addMappingInfo(xml_root, statusOfInterest, fileName)

    # locate intracellular layer in xml tree
    cell_defs = xml_root.find("cell_definitions")
    cell_def = cell_defs.find("cell_definition")
    phenotype = cell_def.find("phenotype")
    intracellular = phenotype.find("intracellular")

    # modify xml bnd_filename, cfg_filename parameters
    bnd_filename = intracellular.find("bnd_filename")
    bnd_filename.text = "config/" + fileName + ".bnd"
    cfg_filename = intracellular.find("cfg_filename")
    cfg_filename.text = "config/" + fileName + ".cfg"

    # modify variable name
    microenvironment = xml_root.find("microenvironment_setup")
    compound = microenvironment.find("variable")
    compound.attrib = {"name": fileName, "units": "dimensionless", "ID": "1"}
    
    # modify decay rate
    #physical_parameter_set = compound.find("physical_parameter_set")
    #decay_rate = physical_parameter_set("decay_rate")
    # MODIFY DECAY RATE HERE

    # modify substrate name (in save)
    save = xml_root.find("save")
    svg = save.find("SVG")
    plot_substrate = svg.find("plot_substrate")
    substrate = plot_substrate.find("substrate")
    substrate.text = fileName

    # modify substrate name (in motility)
    motility = phenotype.find("motility")
    options = motility.find("options")
    chemotaxis = options.find("chemotaxis")
    substrate = chemotaxis.find("substrate")
    substrate.text = fileName
    advanced_chemotaxis = options.find("advanced_chemotaxis")
    chemotactic_sensitivities = advanced_chemotaxis.find("chemotactic_sensitivities")
    chemotactic_sensitivity = chemotactic_sensitivities.find("chemotactic_sensitivity")
    chemotactic_sensitivity.attrib = {"substrate": fileName}

    # modify substrate name (in secretion)
    secretion = phenotype.find("secretion")
    substrate = secretion.find("substrate")
    substrate.attrib = {"name": fileName}

    # modify substrate_name (in user parameters)
    user_params = xml_root.find("user_parameters")
    substrate_name = user_params.find("substrate_name")
    substrate_name.text = fileName

    # add user parameters (dose_interval, drug_amount, add_compound, substrate_name)
    #user_params = xml_root.find("user_parameters") 
    #dose_interval = ET.SubElement(user_params, "dose_interval", {"type": "double", "units": "minutes", "description": "interval between doses of compound"})
    #dose_interval.text = "1440.0"
    #drug_amount = ET.SubElement(user_params, "drug_amount", {"type": "double", "units": "dimensionless", "description": "concentration/amount of compound to add at dosing interval"})
    #drug_amount.text = "50.0"
    #add_compound = ET.SubElement(user_params, "add_compound", {"type": "bool", "units": "", "description": ""})
    #add_compound.text = "true"
    #substrate_name = ET.SubElement(user_params, "substrate_name", {"type": "string", "units":"", "description": ""})
    #substrate_name.text = fileName

    #<dose_interval type="double" units="minutes" description="interval between doses of compound">1440.0</dose_interval>
    #<drug_amount type="double" units="dimensionless" description="concentration/amount of compound to add at dosing interval">50.0</drug_amount>
    #<add_compound type="bool" units="" description="">true</add_compound>

    # create new xml file
    tree.write(fileName + ".xml")
    print("Created file " + fileName + ".xml")

    # create new BND file
    newFileName  = fileName + ".bnd"
    newFile = open(newFileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")
    
    print("Created file " + newFileName)

    # create cfg file
    createCFG(fileName)
    print("Created file " + fileName + ".cfg")

print("--------------------------")
print("Finished with IBMFA files.")
print("--------------------------")

# create files for Stable Motifs results

# change directory --> EDIT TO BE THE LOCATION YOU WANT TO STORE STABLE MOTIFS FILES
os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\Stable Motifs Files")

for i in range(len(stableMotifs)):
    # get base bnd file info
    nodeDict = getBNDdata()

    # load and parse base xml file
    base_xml = open("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\base_model_file.xml")
    tree = ET.parse(base_xml)
    xml_root = tree.getroot()

    # get node(s) of interest and desired state (from csv)
    intervention = (stableMotifs.iloc[i, 0]).replace(" ", "")
    subinterventions = intervention.split("&")
    fileName = ""

    for s in subinterventions:
        nodeOfInterest = s.split("-")[0]
        statusOfInterest = s.split("-")[1]

        # create new node for intervention
        newName = "pro_" + nodeOfInterest
        nodeDict[newName] = "  logic = " + newName  

        # find node logic corresponding to node of interest
        newLogic = nodeDict.get(nodeOfInterest) + ")"

        if(statusOfInterest == "0"):
            # suppression
            newLogic = newLogic[:-1] + " & !" + newName
            physiName = "anti_" + nodeOfInterest
            fileName += physiName
        else:
            # promotion
            newLogic = newLogic[:-1] + " | " + newName
            physiName = "pro_" + nodeOfInterest
            fileName += physiName

        # edit logic for node of interest
        nodeDict[nodeOfInterest] = newLogic

    # add xml mapping data
    addMappingInfo(xml_root, statusOfInterest, fileName)

    # locate intracellular layer in xml tree
    cell_defs = xml_root.find("cell_definitions")
    cell_def = cell_defs.find("cell_definition")
    phenotype = cell_def.find("phenotype")
    intracellular = phenotype.find("intracellular")

    # modify xml bnd_filename, cfg_filename parameters
    bnd_filename = intracellular.find("bnd_filename")
    bnd_filename.text = "config/" + fileName + ".bnd"
    cfg_filename = intracellular.find("cfg_filename")
    cfg_filename.text = "config/" + fileName + ".cfg"

    # modify variable name
    microenvironment = xml_root.find("microenvironment_setup")
    compound = microenvironment.find("variable")
    compound.attrib = {"name": fileName, "units": "dimensionless", "ID": "1"}

    # modify decay rate
    #physical_parameter_set = compound.find("physical_parameter_set")
    #decay_rate = physical_parameter_set("decay_rate")
    # MODIFY DECAY RATE HERE

    # modify substrate name (in motility)
    motility = phenotype.find("motility")
    options = motility.find("options")
    advanced_chemotaxis = options.find("advanced_chemotaxis")
    chemotactic_sensitivities = advanced_chemotaxis.find("chemotactic_sensitivities")
    chemotactic_sensitivity = chemotactic_sensitivities.find("chemotactic_sensitivity")
    chemotactic_sensitivity.attrib = {"substrate": fileName}

    # modify substrate name (in secretion)
    secretion = phenotype.find("secretion")
    substrate = secretion.find("substrate")
    substrate.attrib = {"name": fileName}

    # modify substrate name (in user parameters)
    user_params = xml_root.find("user_parameters")
    substrate_name = user_params.find("substrate_name")
    substrate_name.text = fileName

    # add user parameters (dose_interval, drug_amount, and add_compound)
    #user_params = xml_root.find("user_parameters")
    #dose_interval = ET.SubElement(user_params, "dose_interval", {"type": "double", "units": "minutes", "description": "interval between doses of compound"})
    #dose_interval.text = "1440.0"
    #drug_amount = ET.SubElement(user_params, "drug_amount", {"type": "double", "units": "dimensionless", "description": "concentration/amount of compound to add at dosing interval"})
    #drug_amount.text = "50.0"
    #add_compound = ET.SubElement(user_params, "add_compound", {"type": "bool", "units": "", "description": ""})
    #add_compound.text = "true"
    #substrate_name = ET.SubElement(user_params, "substrate_name", {"type": "string", "units":"", "description": ""})
    #substrate_name.text = fileName

    # create new xml file
    tree.write(fileName + ".xml")
    print("Created file " + fileName + ".xml")

    # create new BND file
    newFileName  = fileName + ".bnd"
    newFile = open(newFileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")
    
    print("Created file " + newFileName)

    # create cfg file
    createCFG(fileName)
    print("Created file " + fileName + ".cfg")

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

    # load and parse base xml file
    base_xml = open("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Variant Model Files\\base_model_file.xml")
    tree = ET.parse(base_xml)
    xml_root = tree.getroot()

    # update logic of nodes for intervention
    # create a new node for pro_source
    newName = "pro_" + source
    nodeDict[newName] = "  logic = " + newName

    # update existing logic for target node
    if (action == "0"):
        # suppression
        nodeDict[target] = nodeDict[target] + " & !" + newName
    else:
        # excitation
        nodeDict[target] = nodeDict[target] + " | " + newName
    
    # add xml mapping data
    physiName = source + "_" + target + "_" + action
    addMappingInfo(xml_root, action, physiName)

    # locate intracellular layer in xml tree
    cell_defs = xml_root.find("cell_definitions")
    cell_def = cell_defs.find("cell_definition")
    phenotype = cell_def.find("phenotype")
    intracellular = phenotype.find("intracellular")

    # modify xml bnd_filename, cfg_filename parameters
    bnd_filename = intracellular.find("bnd_filename")
    bnd_filename.text = "config/" + physiName + ".bnd"
    cfg_filename = intracellular.find("cfg_filename")
    cfg_filename.text = "config/" + physiName + ".cfg"

    # modify variable name
    microenvironment = xml_root.find("microenvironment_setup")
    compound = microenvironment.find("variable")
    compound.attrib = {"name": physiName, "units": "dimensionless", "ID": "1"}

    # modify decay rate
    #physical_parameter_set = compound.find("physical_parameter_set")
    #decay_rate = physical_parameter_set("decay_rate")
    # MODIFY DECAY RATE HERE

    # modify substrate name (in motility)
    motility = phenotype.find("motility")
    options = motility.find("options")
    advanced_chemotaxis = options.find("advanced_chemotaxis")
    chemotactic_sensitivities = advanced_chemotaxis.find("chemotactic_sensitivities")
    chemotactic_sensitivity = chemotactic_sensitivities.find("chemotactic_sensitivity")
    chemotactic_sensitivity.attrib = {"substrate": physiName}

    # modify substrate name (in secretion)
    secretion = phenotype.find("secretion")
    substrate = secretion.find("substrate")
    substrate.attrib = {"name": physiName}

    # modify substrate name (in user parameters)
    user_params = xml_root.find("user_parameters")
    substrate_name = user_params.find("substrate_name")
    substrate_name.text = physiName

    # add user parameters (dose_interval, drug_amount, and add_compound)
    #user_params = xml_root.find("user_parameters")
    #dose_interval = ET.SubElement(user_params, "dose_interval", {"type": "double", "units": "minutes", "description": "interval between doses of compound"})
    #dose_interval.text = "1440.0"
    #drug_amount = ET.SubElement(user_params, "drug_amount", {"type": "double", "units": "dimensionless", "description": "concentration/amount of compound to add at dosing interval"})
    #drug_amount.text = "50.0"
    #add_compound = ET.SubElement(user_params, "add_compound", {"type": "bool", "units": "", "description": ""})
    #add_compound.text = "true"
    #substrate_name = ET.SubElement(user_params, "substrate_name", {"type": "string", "units":"", "description": ""})
    #substrate_name.text = physiName

    # create new xml file
    tree.write(physiName + ".xml")
    print("Created file " + physiName + ".xml")

    # create new BND file
    newFileName = physiName + ".bnd"
    newFile = open(newFileName, "w")
    for n in nodeDict:
        nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n  rate_up = @logic ? $time_scale: 0;\n  rate_down = @logic ? 0 : $time_scale;\n}"
        newFile.write(nodeString + "\n\n")

    print("Created file " + newFileName)

    # create cfg file
    createCFG(physiName)
    print("Created file " + physiName + ".cfg")

print("----------------------------------------------")
print("Finished with Single Edge Perturbations Files.")
print("----------------------------------------------")
