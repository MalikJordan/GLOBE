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
file = 'bfm17.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

light_lim = np.zeros(parameters["simulation"]["iters"])
irrad = np.zeros(parameters["simulation"]["iters"])
zoo_resp = np.zeros(parameters["simulation"]["iters"])
# prod = np.zeros(parameters["simulation"]["iters"])
# exud = np.zeros(parameters["simulation"]["iters"])
fI = np.zeros(parameters["simulation"]["num_days"])
eir = np.zeros(parameters["simulation"]["num_days"])
gpp = np.zeros(parameters["simulation"]["num_days"])
exu = np.zeros(parameters["simulation"]["num_days"])
rspz = np.zeros(parameters["simulation"]["num_days"])
j=0
# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
for iter in range(0,parameters["simulation"]["iters"]-1):
    light_lim[iter], irrad[iter], zoo_resp[iter] = bgc_rate_eqns(iter, base_element, parameters, tracers)
    
    if (iter+1)%(86400/parameters["simulation"]["timestep"]) == 0 and iter > 0:
        fI[j] = np.average(light_lim[240*j:240*(j+1)])
        eir[j] = np.average(irrad[240*j:240*(j+1)])
        gpp[j] = np.average(tracers["phyto1"].gpp[240*j:240*(j+1)])
        exu[j] = np.average(tracers["phyto1"].exu[240*j:240*(j+1)])
        rspz[j] = np.average(zoo_resp[240*j:240*(j+1)])
        j += 1


# x = np.linspace(0,15552000,7)
# marks = ['J','F','M','A','M','J','J']

x = np.linspace(0,62208000,9)
marks = ['J','A','J','O','J','A','J','O','J']

fig, ax = plt.subplots()
ax.plot(np.linspace(0,parameters["simulation"]["num_days"]-1,parameters["simulation"]["num_days"]),fI)
plt.savefig("fI.jpg")

fig, ax = plt.subplots()
ax.plot(np.linspace(0,parameters["simulation"]["num_days"]-1,parameters["simulation"]["num_days"]),eir)
plt.savefig("eir.jpg")

fig, ax = plt.subplots()
ax.plot(np.linspace(0,parameters["simulation"]["num_days"]-1,parameters["simulation"]["num_days"]),gpp)
plt.savefig("gpp.jpg")

fig, ax = plt.subplots()
ax.plot(np.linspace(0,parameters["simulation"]["num_days"]-1,parameters["simulation"]["num_days"]),exu)
plt.savefig("exu.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],light_lim)
plt.title("Light Limitation")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("light_lim.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[3,...])
plt.title("Phyto Chlorophyll-a")
# plt.xlabel("Time [s]")
plt.ylabel("mg Chl a / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("phyto_chl.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[1,...])
# plt.title("Phyto N")
# plt.savefig("phyto_n.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[0,...])
# plt.title("Zoo Carbon")
# # plt.xlabel("Time [s]")
# plt.ylabel("mg C a / m3")
# plt.xlim([0,parameters["simulation"]["time"][-1]])
# # plt.xticks(x,marks)
# plt.savefig("zoo_c.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["dom1"].conc[0,...])
# plt.title("DOC")
# # plt.xlabel("Time [s]")
# plt.ylabel("mg C / m3")
# plt.xlim([0,parameters["simulation"]["time"][-1]])
# # plt.xticks(x,marks)
# plt.savefig("doc.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["pom1"].conc[1,...])
# plt.title("PON")
# # plt.xlabel("Time [s]")
# plt.ylabel("mmol N / m3")
# plt.xlim([0,parameters["simulation"]["time"][-1]])
# # plt.xticks(x,marks)
# plt.savefig("pon.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["o2"].conc[0,...])
plt.title("Oxygen")
# plt.xlabel("Time [s]")
plt.ylabel("mmol O2 / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("oxy.jpg")


fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["no3"].conc[0,...])
plt.title("Nitrate")
# plt.xlabel("Time [s]")
plt.ylabel("mmol N / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("no3.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["po4"].conc[0,...])
plt.title("Phosphate")
# plt.xlabel("Time [s]")
plt.ylabel("mmol P / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("po4.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["nh4"].conc[0,...])
plt.title("Ammonium")
# plt.xlabel("Time [s]")
plt.ylabel("mmol N / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("nh4.jpg")



# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],tracers["no3"].conc[0,...])
# ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[1,...])
# ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[1,...])
# ax.plot(parameters["simulation"]["time"],tracers["pom1"].conc[1,...])
# plt.title("Nitrate")
# plt.xlabel("Time [s]")
# plt.ylabel("mmol N / m3")
# plt.legend(["N","P","Z","D"])
# plt.savefig("nitrogen.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[0,...])
ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[0,...])
plt.title("Carbon")
plt.xlabel("Time [s]")
plt.ylabel("mg C / m3")
ax.legend(["P","Z"])
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("zooc_vs_phytoc.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[1,...])
ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[1,...])
plt.title("Nitrogen")
plt.xlabel("Time [s]")
plt.ylabel("mmol N / m3")
ax.legend(["P","Z"])
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
plt.savefig("zoon_vs_phyton.jpg")

fig, ax = plt.subplots()
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[2,...])
ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[2,...])
plt.title("Phosphorus")
plt.xlabel("Time [s]")
plt.ylabel("mmol P / m3")
plt.xlim([0,parameters["simulation"]["time"][-1]])
# plt.xticks(x,marks)
ax.legend(["P","Z"])
plt.savefig("zoop_vs_phytop.jpg")

# fig, ax = plt.subplots()
# ax.plot(parameters["simulation"]["time"],graze_sum)
# plt.savefig("graze_sum.jpg")