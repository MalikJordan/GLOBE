import os
# import sys
import numpy as np
from setup.initialize import import_bgc_model, import_physical_model
from functions.bgc_rate_eqns import bgc_rate_eqns
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------
# Import and initialize model
# ----------------------------------------------------------------------------------------------------
# check_file = False
# first_check = True
# while not check_file:
#     if first_check:
#         file = input("\nEnter the name of the yaml input file. Press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     else:
#         file = input("\nInput file not found. Re-enter the name of the yaml input file or press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     if file == 'Q' or file == 'q':
#         sys.exit()
#     check_file = os.path.exists(file)
#     first_check = False

# Import physical model
file = 'physical.yaml'
file_path = os.getcwd() + '/' + file
physical = import_physical_model(file_path)

file = 'tests/bfm17/bfm17-1d.yaml'
file_path = os.getcwd() + '/' + file
base_element, reactions, tracers = import_bgc_model(file_path, physical)

# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
# for iter in range(0,physical["simulation"]["iters"]-1):
#     bgc_rate_eqns(iter, base_element, parameters, tracers)


for iter in range(0,physical["simulation"]["iters"]-1):
    t=0












concentration = []
tracer_indices = {}
index = 0   # counting number to keep track of tracer index in concentration matrix
for t in tracers:
    num_constituents = len(tracers[t].composition)
    tracer_indices[t] = list(np.arange(index, index+num_constituents, 1))
    for i in range(num_constituents):
        concentration.append(tracers[t].conc[i,...])

        index += 1

concentration = np.array(concentration,dtype=float)

np.savez('npzd.npz',concentration=concentration,time=physical["simulation"]["time"])
np.savez('tracer_indices_npzd.npz',**tracer_indices)

print('Simulation complete.')
