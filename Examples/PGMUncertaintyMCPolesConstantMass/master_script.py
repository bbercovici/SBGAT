import os
import json
import numpy as np

import os
import platform
import sys
import itertools
import time
import socket 

def generate_all_cases_dictionnary_list(base_dictionnary,all_cases_dictionnary,base_location,sim_name):

    if len(all_cases_dictionnary) > 0:
        keys, values = zip(*all_cases_dictionnary.items())
        dictionnary_list = [dict(zip(keys, v)) for v in itertools.product(*values)]
        all_cases_dictionnary_list = [{**dictionnary_list[e],**base_dictionnary} for e in range(len(dictionnary_list))]
    else:
        all_cases_dictionnary_list = [base_dictionnary]
        dictionnary_list = [base_dictionnary]
    for e in range(len(dictionnary_list)):
        all_cases_dictionnary_list[e]["INPUT_DIR"] = base_location + "input/" + sim_name + "_" + str(e) + "/"
        all_cases_dictionnary_list[e]["OUTPUT_DIR"] = base_location + "output/" + sim_name + "_" + str(e) + "/"

    return all_cases_dictionnary_list


# Replace the paths after 'base_location' with the existing directory under which the input/ and /output sub-directories
# will be created and populated
if (socket.gethostname() == "fortuna"):
    base_location = "/orc_raid/bebe0705/PGMUncertaintyMCPolesConstantMass/"
else:
    base_location = "../"

# SIM_PREFIX will be added to the name of every folder to be put in input/ and output/ 
SIM_PREFIX = "PGMUncertaintyMCPolesConstantMass"

# Dictionnary storing simulation inputs to be kept constant
base_dictionnary = {
"DENSITY" : 4500,
"UNIT_IN_METERS" : False,
"STEP_SIZE" : 10e3,
"PATH_SHAPE" : "../../../resources/shape_models/psyche.obj",
"PROJECTION_AXIS" : 0,
"UNCERTAINTY_TYPE" : "normal",
"N_MONTE_CARLO" : 5000,
"CORRELATION_DISTANCE" : 75e3,
"ERROR_STANDARD_DEV" : 10e3,
"COV_REGION_CENTERS" : [0,1147]
}

# Dictionnary storing simulation inputs to be looped over
# for instance, one could have
# all_cases_dictionnary = {
# "INPUT_THAT_MUST_BE_CHANGED_1" : [1,2,3]
# "INPUT_THAT_MUST_BE_CHANGED_2" : [True,False]
# }
# which means that a total of six (two times three) simulations will be run
# and saved in input/ and output/, with the names of the subfolder prefixed by SIM_PREFIX"

all_cases_dictionnary = {
"HOLD_MASS_CONSTANT" : [True]
}

# There shouldn't be any reason to modify the following
all_data = generate_all_cases_dictionnary_list(base_dictionnary,
all_cases_dictionnary,base_location,SIM_PREFIX)
os.system("cmake .. && make")
for data in all_data:

    print("\t Case " + data["INPUT_DIR"].split("/")[-2])
    print("\t - Making directory")

    os.system("mkdir " + data["INPUT_DIR"])
    os.system("mkdir " + data["OUTPUT_DIR"])

    print("\t - Copying input file in build/")

    with open('input_file.json', 'w') as outfile:
        json.dump(data, outfile)
    print("\t - Saving input file in input/ and output/")
    with open(data["INPUT_DIR"] + 'input_file.json', 'w') as outfile:
        json.dump(data, outfile)
    with open(data["OUTPUT_DIR"] + 'input_file.json', 'w') as outfile:
        json.dump(data, outfile)
    print("\t - Running case " +  data["INPUT_DIR"].split("/")[-2])
    os.system("> " + data["OUTPUT_DIR"] + "log.txt")
    os.system("./PGMUncertaintyMCPolesConstantMass 2>&1 | tee " + data["OUTPUT_DIR"] + "log.txt" )
   

