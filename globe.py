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
plt.plot(time,solution[1,:])
plt.plot(time,solution[2,:])
plt.plot(time,solution[3,:])
plt.plot(time,solution[0,:])
# plt.legend(['DIN', 'Phytoplankton'], loc='upper left')
# plt.legend(['DIN', 'Phytoplankton', 'Zooplankton'], loc='upper left')
plt.legend(['DIN', 'Phytoplankton', 'Zooplankton', 'Detritus'], loc='upper left')
ticks = [15,105,195,285,375,465,555,645,735,825,915,1005]
for i in range(0,len(ticks)):
    ticks[i] = ticks[i] * 86400

plt.xticks(ticks, ['J','A','J','S','J','A','J','S','J','A','J','S'])
plt.xlim([0,end_time])
# plt.ylim([0,10.5])
plt.ylabel("mmol N / m^3")
plt.savefig('npzd-iterative.jpg')

print()