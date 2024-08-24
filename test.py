from functions.seasonal_cycling import get_temperature, get_mixed_layer_depth, get_sunlight, get_salinity
from setup.initialize import coordinate_system
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def test_coordinate_system():

    # Open json file containing reactions and model parameters
    file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        parameters = model_info["parameters"]

    # Extract model parameters
    water_column_parameters = parameters["water_column"]

    water_column_parameters = coordinate_system(water_column_parameters)

    return water_column_parameters
    


def test_get_temperature(time):
    temp = np.zeros_like(time)

    # Open json file containing reactions and model parameters
    file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        parameters = model_info["parameters"]

    # Extract model parameters
    environmental_parameters = parameters["environment"]

    for i in range(0,len(time)):
        temp[i] = get_temperature(time[i], environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])

    fig1, ax1 = plt.subplots()
    ax1.plot(time, temp)
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Temperature [C]')
    plt.savefig('test_get_temperature.jpg')

def test_get_mld(time):
    mld = np.zeros_like(time)

    # Open json file containing reactions and model parameters
    file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        parameters = model_info["parameters"]

    # Extract model parameters
    environmental_parameters = parameters["environment"]

    for i in range(0,len(time)):
        mld[i] = get_mixed_layer_depth(time[i], environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])

    fig2, ax2 = plt.subplots()
    ax2.plot(time, mld)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Mixed Layer Depth [m]')
    plt.savefig('test_get_mld.jpg')


def test_get_sunlight(time):
    sun = np.zeros_like(time)

    # Open json file containing reactions and model parameters
    file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        parameters = model_info["parameters"]

    # Extract model parameters
    environmental_parameters = parameters["environment"]

    for i in range(0,len(time)):
        sun[i] = get_sunlight(time[i], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])

    fig3, ax3 = plt.subplots()
    ax3.plot(time, sun)
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('Surface Irradiance')
    plt.savefig('test_get_sunlight.jpg')


def tets_get_salinity(time):
    salt = np.zeros_like(time)

    # Open json file containing reactions and model parameters
    file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        parameters = model_info["parameters"]

    # Extract model parameters
    environmental_parameters = parameters["environment"]

    for i in range(0,len(time)):
        salt[i] = get_salinity(time[i], environmental_parameters["winter_salt"], environmental_parameters["summer_salt"])

    fig3, ax3 = plt.subplots()
    ax3.plot(time, salt)
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('Salinity')
    plt.savefig('test_get_salinity.jpg')


dt = 360
seconds_per_day = 86400
simulation_length = 30     # days
end = seconds_per_day * simulation_length
steps = end/dt + 1
time = np.linspace(0,int(end),int(steps))

# test_get_temperature(time)
# test_get_sunlight(time)

dt = 360
seconds_per_day = 86400
simulation_length = 360     # days
end = seconds_per_day * simulation_length
steps = end/dt + 1
time = np.linspace(0,int(end),int(steps))
# test_get_mld(time)
# tets_get_salinity(time)

# water_column_parameters = test_coordinate_system()
print()


max_growth_rate = 1. / seconds_per_day              # s^-1
nitrogen_half_saturation = 1.                       # μmol N l^-1
max_grazing_rate = 1. / seconds_per_day             # s^-1
zooplankton_death_rate = 0.2 / seconds_per_day      # s^-1
grazing = 0.2                                       # μmol N l^-1
phytoplankton_death_rate = 0.1 / seconds_per_day    # s^-1
zooplankton_assimiilated_nitrogen = 0.7
light_intensity = 0.25

iterations_needed = simulation_length * seconds_per_day / dt
t_span = np.linspace(0,int(iterations_needed),int(iterations_needed)+1)

def dYdt(Y,t):
    p, z, n = Y
    dpdt = ( max_growth_rate * (n/ (nitrogen_half_saturation + n)) * light_intensity * p ) \
                - ( max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) - ( phytoplankton_death_rate * p )
    dzdt = ( zooplankton_assimiilated_nitrogen * max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) - ( zooplankton_death_rate * z )
    dndt = -( max_growth_rate * (n/ (nitrogen_half_saturation + n)) * light_intensity * p ) \
                - ( (1-zooplankton_assimiilated_nitrogen) * max_grazing_rate * (1 - np.exp(-grazing*p)) * z ) \
                    + ( phytoplankton_death_rate * p ) + ( zooplankton_death_rate * z )

    return [dpdt, dzdt,  dndt]

# t_span = np.linspace(0,int(params_POMBFM.idays)*seconds_per_day)
p0 = 2.5
z0 = 1.5
n0 = 4.
Y0 = [p0, z0, n0]

sol = odeint(dYdt, Y0, t_span)
# sol2 = solve_ivp(dYdt, t_span, Y0)

fig1, ax = plt.subplots()
# fig1 = plt.figure()
# ax = fig1.add_subplot(1,3,1)
phyto, = ax.plot(t_span,sol[:,0], label='Phytoplankton')
zoo, = ax.plot(t_span,sol[:,1], label='Zooplankton')
nut, = ax.plot(t_span,sol[:,2], label='Nutrients')
ax.legend(loc='upper right')
plt.title('odeint')
plt.xlabel('Time (seconds)')
plt.ylabel('Concentration (μmol N l^-1)')
plt.savefig('test_npz.jpg')

# fig1, ax = plt.subplots()
# ax = fig1.add_subplot(1,3,2)
# phyto, = ax.plot(t_span,sol2[0,:], label='Phytoplankton')
# zoo, = ax.plot(t_span,sol2[1,:], label='Zooplankton')
# nut, = ax.plot(t_span,sol2[2,:], label='Nutrients')
# ax.legend(loc='upper right')
# plt.title('solve_ivp')
# plt.xlabel('Time (seconds)')
# plt.ylabel('Concentration (μmol Ν l^-1)')
# plt.show()

# ---------- 0D ----------
NPZ = np.zeros((3,int(iterations_needed)+1))
NPZ[0,0] = 2.5    # μmol N l^-1
NPZ[1,0] = 1.5    # μmol N l^-1
NPZ[2,0] = 4.     # μmol N l^-1

for i in range(0,int(iterations_needed)):
    phyto_coeff = max_growth_rate * ( NPZ[2,i] / (nitrogen_half_saturation + NPZ[2,i]) ) * light_intensity
    zoo_coeff1 = max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    zoo_coeff2 = zooplankton_assimiilated_nitrogen * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    zoo_coeff3 = ( 1 - zooplankton_assimiilated_nitrogen ) * max_grazing_rate * ( 1 - np.exp(-grazing*NPZ[0,i]) )
    
    NPZ[0,i+1] = NPZ[0,i] + dt * ( phyto_coeff * NPZ[0,i] - zoo_coeff1 * NPZ[1,i] - phytoplankton_death_rate * NPZ[0,i] )
    NPZ[1,i+1] = NPZ[1,i] + dt * ( zoo_coeff2 * NPZ[1,i] - zooplankton_death_rate * NPZ[1,i] )
    NPZ[2,i+1] = NPZ[2,i] + dt * ( -phyto_coeff * NPZ[0,i] + zoo_coeff3 * NPZ[1,i] + phytoplankton_death_rate * NPZ[0,i] + zooplankton_death_rate * NPZ[1,i] )
    
fig2, ax = plt.subplots()
# ax = fig1.add_subplot(1,3,3)
phyto, = ax.plot(t_span,NPZ[0,:], label='Phytoplankton')
zoo, = ax.plot(t_span,NPZ[1,:], label='Zooplankton')
nut, = ax.plot(t_span,NPZ[2,:], label='Nutrients')
ax.legend(loc='upper right')
plt.title('forward euler')
plt.xlabel('Time (days)')
plt.ylabel('Concentration (μmol N l^-1)')
plt.savefig('test_npz-iterative.jpg')

