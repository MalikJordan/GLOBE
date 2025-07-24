import os
# import sys
import numpy as np
from setup.initialize import import_model
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
# file = 'bfm17.yaml'
file = 'bfm17-OOI.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

j=0
irrad = np.zeros(parameters["simulation"]["iters"])
light_lim = np.zeros(parameters["simulation"]["iters"])
nut_lim = np.zeros(parameters["simulation"]["iters"])
# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
for iter in range(0,parameters["simulation"]["iters"]-1):
    # light_lim[iter], irrad[iter], nut_lim[iter] = bgc_rate_eqns(iter, base_element, parameters, tracers)
    bgc_rate_eqns(iter, base_element, parameters, tracers)
    # bgc_rate_eqns(iter)

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

np.savez('bfm17-OOI.npz',concentration=concentration,time=parameters["simulation"]["time"])
np.savez('tracer_indices_bfm17-OOI.npz',**tracer_indices)

print('Simulation complete.')

# np.savez('phyto-0424.npz',gpp=tracers['phyto1'].gpp, rsp=tracers['phyto1'].rsp, exu=tracers['phyto1'].exu, 
#          lys=tracers['phyto1'].lys, npp=tracers['phyto1'].npp, psn=tracers['phyto1'].psn,
#          uptn=tracers['phyto1'].uptn, uptp=tracers['phyto1'].uptp)

# np.savez('light_globe.npz',irrad=irrad, light_lim=light_lim)
# np.savez('nut_lim_globe.npz',nut_lim=nut_lim)
x=1

# # Convert time array from seconds to months
# sec_mon = 60 * 60 * 24 * 30
# time = parameters["simulation"]["time"]/sec_mon

# # Ticks and labels for plots
# xticks = [0,2,4,6,8,10,12]
# xlabel = ['J','M','M','J','S','N','J']

# # Create plots

# # Folder path
# folder = os.getcwd() + '/tests/2npzd-0404'

# # Nutrients
# fig, axs = plt.subplots(1,2,figsize=(12,5))

# axs[0].plot(time,concentration[0])
# axs[0].set_title("Nitrate")
# axs[0].set_xlabel("Time [months]")
# axs[0].set_xticks(xticks,xlabel)
# axs[0].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1].plot(time,concentration[1])
# axs[1].set_title("Phosphate")
# axs[1].set_xlabel("Time [months]")
# axs[1].set_xticks(xticks,xlabel)
# axs[1].set_ylabel("mmol P $\mathregular{m^{-3}}$")

# fig.suptitle("Nutrients")
# fig.tight_layout()
# nut = os.path.join(folder,"nutrients.jpg")
# plt.savefig(nut)

# # Phytoplankton
# fig, axs = plt.subplots(1,2,figsize=(12,5))

# axs[0].plot(time,concentration[2])
# axs[0].set_title("Nitrate")
# axs[0].set_xlabel("Time [months]")
# axs[0].set_xticks(xticks,xlabel)
# axs[0].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1].plot(time,concentration[3])
# axs[1].set_title("Phosphate")
# axs[1].set_xlabel("Time [months]")
# axs[1].set_xticks(xticks,xlabel)
# axs[1].set_ylabel("mmol P $\mathregular{m^{-3}}$")

# fig.suptitle("Phytoplankton")
# fig.tight_layout()
# phyto = os.path.join(folder,"phytoplankton.jpg")
# plt.savefig(phyto)

# # Zooplankton
# fig, axs = plt.subplots(1,2,figsize=(12,5))

# axs[0].plot(time,concentration[4])
# axs[0].set_title("Nitrate")
# axs[0].set_xlabel("Time [months]")
# axs[0].set_xticks(xticks,xlabel)
# axs[0].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1].plot(time,concentration[5])
# axs[1].set_title("Phosphate")
# axs[1].set_xlabel("Time [months]")
# axs[1].set_xticks(xticks,xlabel)
# axs[1].set_ylabel("mmol P $\mathregular{m^{-3}}$")

# fig.suptitle("Zooplankton")
# fig.tight_layout()
# zoo = os.path.join(folder,"zooplankton.jpg")
# plt.savefig(zoo)

# # Detritus
# fig, axs = plt.subplots(1,2,figsize=(12,5))

# axs[0].plot(time,concentration[6])
# axs[0].set_title("Nitrate")
# axs[0].set_xlabel("Time [months]")
# axs[0].set_xticks(xticks,xlabel)
# axs[0].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1].plot(time,concentration[7])
# axs[1].set_title("Phosphate")
# axs[1].set_xlabel("Time [months]")
# axs[1].set_xticks(xticks,xlabel)
# axs[1].set_ylabel("mmol P $\mathregular{m^{-3}}$")

# fig.suptitle("Detritus")
# fig.tight_layout()
# pom = os.path.join(folder,"detritus.jpg")
# plt.savefig(pom)
