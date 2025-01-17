import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# SETUP_PATH = Path()
from setup.initialize import import_model
from functions.bgc_rate_eqns import bgc_rate_eqns

# ----------------------------------------------------------------------------------------------------
# GLOBE model simulation
# ----------------------------------------------------------------------------------------------------
# Import and initialize model
file = 'test_npz.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

# Begin simulation
for iter in range(0,parameters["simulation"]["iters"]-1):
    bgc_rate_eqns(iter, base_element, parameters, tracers)

fig, ax = plt.subplots()
x = [0,2592000,5184000,7776000,10268000,12960000,15552000]
marks = ['J','F','M','A','M','J','J']
ax.plot(parameters["simulation"]["time"],tracers["no3"].conc[0,:])
ax.plot(parameters["simulation"]["time"],tracers["phyto1"].conc[0,:])
ax.plot(parameters["simulation"]["time"],tracers["zoo1"].conc[0,:])
plt.ylabel("mmol N / m^3")
plt.xlabel('Time [months]')
plt.xlim([0,15552000])
plt.xticks(x,marks)

# Shrink current axis by 10%
box = ax.get_position()
ax.set_position([box.x0+0.04, box.y0+0.02, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(['N','P','Z'],loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("test_npz.jpg")