import numpy as np
import matplotlib.pyplot as plt

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
plt.figure()
plt.plot(x,N)
plt.plot(x,P)
plt.plot(x,Z)
plt.plot(x,D)
plt.legend(['N','P','Z','D'])

plt.savefig('riley.jpg')
