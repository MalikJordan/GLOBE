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
file = '2npzd.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

# Begin simulation
for iter in range(0,parameters["simulation"]["iters"]-1):
    bgc_rate_eqns(iter, base_element, parameters, tracers)

# ----------------------------------------------------------------------------------------------------
# Test case and parameters from Riley Brady
# ----------------------------------------------------------------------------------------------------
# Here we set up the default parameters/coefficients. 
DT = 1 # Time Step (in days)
NUM_STEPS = 365 # Number of time steps to be computed and plotted

# Temperature-Dependent Growth Rate
a = 0.6
b = 1.066
c = 1
T = 15
Vm = a * b**(c*T) # Maximum growth rate (per day)

# Other parameters
Kn = 1    # Half-saturation constant for nitrogen uptake (umolN per l)
Rm = 1    # Maximum grazing rate (per day)
g = 0.2  # Zooplankton death rate (per day)
lambda_Z = 0.2  # Grazing constant (umolN per l)
epsilon = 0.1  # Phyto death rate (per day)
f = 0.25 # Light intensity (assumed constant)

# Detritus-related stuff.
alpha = 0.0 # Fraction of zoo. uptake that goes immediately to dissolved nutrients.
beta = 1.0  # Assimilation efficiency of zooplankton.
r = 0.15 # Respiration rate.
phi = 0.4 # Remineralization rate of detritus.

# Set Initial Conditions (umol per L)
N_0 = 4 
P_0 = 2.5 
Z_0 = 1.5
D_0 = 0

# Initialize Arrays
N = np.empty(NUM_STEPS, dtype="float")
P = np.empty(NUM_STEPS, dtype="float")
Z = np.empty(NUM_STEPS, dtype="float")
D = np.empty(NUM_STEPS, dtype="float")

# Insert Initial Values
N[0] = N_0
P[0] = P_0
Z[0] = Z_0
D[0] = D_0

# Here we use the Euler forward method to solve for t+1 and reference t. 
for idx in np.arange(1, NUM_STEPS, 1):
    t = idx - 1
    
    # Common terms for simpler code
    gamma_N   = N[t] / (Kn + N[t])
    zoo_graze = Rm * (1 - np.exp(-lambda_Z * P[t])) * Z[t]
    
    # Equation calculations
    N[idx] = DT * (-Vm*gamma_N*f*P[t] + alpha*zoo_graze + epsilon*P[t] + g*Z[t] + phi*D[t]) + N[t] 
    P[idx] = DT * (Vm*gamma_N*f*P[t] - zoo_graze - epsilon*P[t] - r*P[t]) + P[t]
    Z[idx] = DT * (beta*zoo_graze - g*Z[t]) + Z[t]  
    D[idx] = DT * (r*P[t] + (1-alpha-beta)*zoo_graze - phi*D[t]) + D[t]

x = np.arange(1, NUM_STEPS + 1, 1)

# ----------------------------------------------------------------------------------------------------
# Plot results
# ----------------------------------------------------------------------------------------------------
fig, ax = plt.subplots()
x = parameters["simulation"]["time"][0:-1]
ax.plot(x,tracers["no3"].conc[0,:-1])
# ax.plot(x,N,'--k')
ax.plot(x,tracers["phyto1"].conc[0,:-1])
# ax.plot(x,P,'--k')
ax.plot(x,tracers["zoo1"].conc[0,:-1])
# ax.plot(x,Z,'--k')
ax.plot(x,tracers["pom1"].conc[0,:-1])
# ax.plot(x,D,'--k')

plt.ylabel("mmol N / m^3")
plt.xlabel("Time (days)")
# plt.xlim([0,365])

# Shrink current axis by 10%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.legend(['N','','P','','Z','','D',''],loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(['N','P','Z','D'],loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("2npzd_n.jpg")


fig, ax = plt.subplots()

ax.plot(x,tracers["po4"].conc[0,:-1])
ax.plot(x,tracers["phyto1"].conc[1,:-1])
ax.plot(x,tracers["zoo1"].conc[1,:-1])
ax.plot(x,tracers["pom1"].conc[1,:-1])

plt.ylabel("mmol P / m^3")
plt.xlabel("Time (days)")
# plt.xlim([0,365])

# Shrink current axis by 10%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(['N','P','Z','D'],loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("2npzd_p.jpg")



