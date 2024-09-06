import os
import numpy as np
from setup.initialize import import_model
from functions.other_functions import *
from functions import rates
from functions.bgc_rate_eqns import bgc_rate_eqns
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------
# Import and initialize model
conc, parameters, tracers = import_model()

# ----------------------------------------------------------------------------------------------------
# Begin simulation
end_time = 86400 * parameters["simulation"]["num_days"]
t_span = (0,end_time)
t_eval = np.arange(0, end_time, parameters["simulation"]["timestep"])

# NP
# conc[0] = 4.0   # DIN
# conc[1] = 2.5   # Phyto

# NPZ
# conc[0] = 4.0   # DIN
# conc[1] = 2.5   # Phyto
# conc[2] = 1.5   # Zoo

# NPZD
# conc[0] = 0.5   # Detritus
# conc[1] = 4.0   # DIN
# conc[2] = 2.5   # Phyto
# conc[3] = 1.5   # Zoo
conc[0] = 0.1   # Detritus
conc[1] = 10.0   # DIN
conc[2] = 0.4   # Phyto
conc[3] = 0.4   # Zoo

# model_solution = solve_ivp(bgc_rate_eqns, t_span, conc, args=(tracers,), method='RK23', t_eval=t_eval)

# ----------------------------------------------------------------------------------------------------
# Plot solution
# fig1, ax1 = plt.subplots()

# NP
# ax1.plot(model_solution.t, model_solution.y[0])
# ax1.plot(model_solution.t, model_solution.y[1])

# ax1.legend(['DIN', 'Phytoplankton'])
# plt.savefig('np.jpg')

# NPZ
# ax1.plot(model_solution.t, model_solution.y[0])
# ax1.plot(model_solution.t, model_solution.y[1])
# ax1.plot(model_solution.t, model_solution.y[2])

# ax1.legend(['DIN', 'Phytoplankton', 'Zooplankton'])
# plt.savefig('npz.jpg')

# NPZD
# ax1.plot(model_solution.t, model_solution.y[1])
# ax1.plot(model_solution.t, model_solution.y[2])
# ax1.plot(model_solution.t, model_solution.y[3])
# ax1.plot(model_solution.t, model_solution.y[0])

# ax1.legend(['DIN', 'Phytoplankton', 'Zooplankton', 'Detritus'])
# plt.savefig('npzd.jpg')

# ----------------------------------------------------------------------------------------------------
# Here we set up the default parameters/coefficients. 
DT = 1 # Time Step (in days)
NUM_STEPS = 180 # Number of time steps to be computed and plotted

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

end_time = 86400 * parameters["simulation"]["num_days"]
time = np.arange(0, end_time + parameters["simulation"]["timestep"], parameters["simulation"]["timestep"])

solution = np.zeros((len(conc),len(time)))
solution[:,0] = conc
for i in range(0,len(time)):
    bgc_rates = bgc_rate_eqns(time[i], conc, tracers)

    conc += bgc_rates * parameters["simulation"]["timestep"]
    solution[:,i] = conc

plt.subplots()
# plt.plot(time,solution[0,:])
plt.plot(time,solution[1,:],'-b')
plt.plot(x,N,'--b')
plt.plot(time,solution[2,:],'y')
plt.plot(x,P,'--y')
plt.plot(time,solution[3,:],'-g')
plt.plot(x,Z,'--g')
plt.plot(time,solution[0,:],'-r')
plt.plot(x,D,'--r')

# plt.legend(['DIN', 'Phytoplankton'], loc='upper left')
# plt.legend(['DIN', 'Phytoplankton', 'Zooplankton'], loc='upper left')
# plt.legend(['DIN', 'Phytoplankton', 'Zooplankton', 'Detritus'], loc='upper left')
# ticks = [15,105,195,285,375,465,555,645,735,825,915,1005]
# for i in range(0,len(ticks)):
#     ticks[i] = ticks[i] * 86400

# plt.xticks(ticks, ['J','A','J','S','J','A','J','S','J','A','J','S'])
plt.xlim([0,end_time])
# plt.ylim([0,10.5])
plt.ylabel("mmol N / m^3")
plt.savefig('npzd-check.jpg')

print()