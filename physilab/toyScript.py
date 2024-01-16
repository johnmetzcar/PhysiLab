import os
import pandas

# set working directory
# os.chdir("C:\\Users\\pletz\\OneDrive\\Desktop\\IU\\Y390\\Script for Toy Model")

# import data from csv
results = pandas.read_csv("toyResults.csv", header=None) # how will we handle more than one intervention per line?

# import the original bnd file
og = open("model_0.bnd", "r+")
bndText = og.readlines()

# create one big long string of the entire bnd file
bndString = ""
for i in range(len(bndText)):
    bndString = bndString + bndText[i]

# separate the string into a list of strings: one for each node
nodes = bndString.split("Node")

# remove the empty string at the beginning
nodes.pop(0)

# store in a better format: make a dictionary with key = node and value = logic 
nodeDict = {}
for node in nodes:
    nodename = node[1:2]
    nodelogic = node.split("\n")[1]
    nodeDict[nodename] = nodelogic

# loop through the interventions
for i in range(len(results)):
    # get node of interest and desired state
    intervention = results.iloc[i][0]
    nodeOfInterest = intervention.split("-")[0]
    #print("Node of Interest: ", nodeOfInterest)
    statusOfInterest = intervention.split("-")[1]
    #print("State of Node of Interest: ", statusOfInterest)

    # create new node for intervention
    newName = "pro_" + nodeOfInterest
    nodeDict[newName] = " logic = " + newName

    # find node logic corresponding to node of interest
    newlogic = nodeDict.get(nodeOfInterest)

    # edit logic
    if statusOfInterest == "1":
      newlogic = newlogic[:-1] + " | " + newName
      nodeDict[nodeOfInterest] = newlogic
    # WILL NEED TO ADD LOGIC FOR SUPPRESSION OF NODE OF INTEREST

# create new BND file
newFile = open("config/updatedBND.bnd", "w") # need to change the output based on the interventions
for n in nodeDict:
    nodeString = "Node " + n + " {\n" + nodeDict.get(n) + ";\n rate_up = @logic ? $time_scale : 0;\n rate_down = @logic ? 0 : $time_scale;\n}"
    newFile.write(nodeString + '\n\n')

print("Done!")