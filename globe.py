import os
import sys
import numpy as np
from setup.initialize import import_model
from functions.other_functions import *
from functions import rates
from functions.bgc_rate_eqns import bgc_rate_eqns
from scipy.integrate import odeint, solve_ivp
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
file = 'test_npzd.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
for iter in range(0,parameters["simulation"]["iters"]-1):
    # bgc_rate_eqns(parameters["simulation"]["time"][i-1], file_path, tracers)
    bgc_rate_eqns(iter, base_element, parameters, tracers)

#     conc += bgc_rates * parameters["simulation"]["timestep"]
#     solution[:,i] = conc


# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[3,...])
# plt.title("Phyto Chlorophyll-a")
# plt.xlabel("Time [s]")
# plt.ylabel("mg Chl a / m3")
# plt.savefig("phyto_chl.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[0,...])
plt.title("Phyto N")
plt.savefig("npzd_phyto.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[0,...])
# plt.title("Zoo Carbon")
# plt.xlabel("Time [s]")
# plt.ylabel("mg C a / m3")
# plt.savefig("zoo_c.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["dom1"].conc[0,...])
# plt.title("DOC")
# plt.xlabel("Time [s]")
# plt.ylabel("mg C / m3")
# plt.savefig("doc.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["o2"].conc[0,...])
# plt.title("Oxygen")
# plt.xlabel("Time [s]")
# plt.ylabel("mmol O2 / m3")
# plt.savefig("oxy.jpg")


fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["no3"].conc[0,...])
plt.title("Nitrate")
plt.xlabel("Time [s]")
plt.ylabel("mmol N / m3")
plt.savefig("no3.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["nh4"].conc[0,...])
# plt.title("Ammonium")
# plt.xlabel("Time [s]")
# plt.ylabel("mmol N / m3")
# plt.savefig("nh4.jpg")
