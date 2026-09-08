import os
import numpy as np
import matplotlib.pyplot as plt

folder = os.getcwd() + '/tests/npzd'

# ----------------------------------------------------------------------------------------------------
# GLOBE model simulation
# ----------------------------------------------------------------------------------------------------
# Get GLOBE Data --------------------------------------------------------------------
# Load solution
# path = os.getcwd() + "/tests/npzd/data/npzd.npz"
path = os.getcwd() + "/concentration_npzd_0907.npz"
solution = np.load(path, allow_pickle=True)
# conc = solution["concentration"]     # concentratrion matrix
conc = solution["daily"]
# time = solution["time"]     # time array

# Load tracer indices
# path = os.getcwd() + "/tests/npzd/data/tracer_indices_npzd.npz"
path = os.getcwd() + "/tests/npzd/data/tracer_indices_npzd.npz"
indices = np.load(path, allow_pickle=True)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

no3 = conc[tracer_indices["no3"][0]][0]
phyto = conc[tracer_indices["phyto1"][0]][0]
zoo = conc[tracer_indices["zoo1"][0]][0]
pom = conc[tracer_indices["pom1"][0]][0]

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

    # bgc_rate_eqns(t, base_element, parameters, tracers)

    # dn = N[idx] - tracers["no3"].conc[0][idx]
    # dp = P[idx] - tracers["phyto1"].conc[0][idx]
    # dz = Z[idx] - tracers['zoo1'].conc[0][idx]
    # dd = D[idx] - tracers['pom1'].conc[0][idx]

    # pause = 1


# x = np.arange(1, NUM_STEPS + 1, 1)
x = np.arange(0, NUM_STEPS, 1)

# ----------------------------------------------------------------------------------------------------
# Plot results
# ----------------------------------------------------------------------------------------------------
# fig, ax = plt.subplots()

months = [0,30,60,90,120,150,180]
marks = ['J','F','M','A','M','J','J']
marks_blank = ['','','','','','','']

fig, axs = plt.subplots(2,2,figsize=(10,10), sharex=True)

# axs[0,0].plot(x,tracers["no3"].conc[0,:-1],label='GLOBE')
axs[0,0].plot(x,no3,label='GLOBE')
axs[0,0].plot(x,N,linestyle=(0, (5, 10)),color='black',label='NPZD')
axs[0,0].set_title("(a) Nitrate")
# axs[0,0].set_xticks(months,marks_blank)
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mmol N ${m^{-3}}$")

# axs[0,1].plot(x,tracers["pom1"].conc[0,:-1],label='GLOBE')
axs[0,1].plot(x,pom,label='GLOBE')
axs[0,1].plot(x,D,linestyle=(0, (5, 10)),color='black',label='NPZD')
axs[0,1].set_title("(b) Particulate Organic Nitrogen")
# axs[0,1].set_xticks(months,marks_blank)
axs[0,1].set_xlim([0,360])

# axs[1,0].plot(x,tracers["phyto1"].conc[0,:-1],label='GLOBE')
axs[1,0].plot(x,phyto,label='GLOBE')
axs[1,0].plot(x,P,linestyle=(0, (5, 10)),color='black',label='NPZD')
axs[1,0].set_title("(c) Phytoplankton")
axs[1,0].set_xlabel("Time [days]")
# axs[1,0].set_xticks(months,marks)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mmol N ${m^{-3}}$")

# axs[1,1].plot(x,tracers["zoo1"].conc[0,:-1],label='GLOBE')
axs[1,1].plot(x,zoo,label='GLOBE')
axs[1,1].plot(x,Z,linestyle=(0, (5, 10)),color='black',label='NPZD')
axs[1,1].set_title("(d) Zooplankton")
axs[1,1].set_xlabel("Time [days]")
# axs[1,1].set_xticks(months,marks)
axs[1,1].set_xlim([0,360])

handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles,labels, loc='lower center', ncol=2)

fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])

npzd = os.path.join(folder + "/figures", "npzd_0907.jpg")
plt.savefig(npzd)

error = np.zeros((4,len(N)))
error[0,:] = (N - no3) / (N + 1.E-20)
error[1,:] = (P - phyto) / (P + 1.E-20)
error[2,:] = (Z - zoo) / (Z + 1.E-20)
error[3,:] = (D - pom) / (D + 1.E-20)
print(np.max(error))
x=1