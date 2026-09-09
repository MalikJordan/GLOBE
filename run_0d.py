import os
import time
import numpy as np
import yaml
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from setup.initialize import import_bgc_model, import_physical_model
from functions.seasonal_cycling import get_mixed_layer_depth, get_salinity, get_sunlight, get_temperature, get_wind
from functions.bgc_rate_eqns import bgc_rate_eqns
from functions.calculate_averages import average
from functions.other_functions import concentration_ratio, light_attenuation
from pom.calculations import density_profile
np.set_printoptions(precision=20)

def create_function_inputs(iters, tracers):
    """
    Definition: Takes tracer dictionary and creates lists, arrays, or typed.Dicts for numba calculations

    :return: concentration (array), sinking velocities (array), tracer map (typed.Dict), tracer types (list)
    """

    # Create list of concentrations
    initial_concentration = []

    # Create typed.Dict of tracer indices in concentration
    tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
    
    # Create list of trcaer types
    tracer_type = []   # used in vertical diffusivity calculations

    # Create list of sinking velocities for each tracer
    sinking = []

    index = 0   # counting number for tracer indices
    for trac in tracers:
        num_constituents = len(tracers[trac].composition)   # number of constituents in tracer

        lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
        for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
        tracer_map[trac] = lst  # identify tracer constituents with their own index

        for i in range(num_constituents):
            # add concentration to matrix
            initial_concentration.append(tracers[trac].initial_conc[i,...])    # add concentration to matrix

            # add tracer type to list
            if tracers[trac].type == "detritus":    tracer_type.append(tracers[trac].form)     # need to distinguish particulate/dissolved form
            else:   tracer_type.append(tracers[trac].type)     # just the type

            # add sinking velocity to list
            if hasattr(tracers[trac],"sinking_velocity"):   sinking.append(tracers[trac].sinking_velocity)
            else:   sinking.append(np.zeros(tracers[trac].initial_conc.shape[1]))

            # add tracer type to list
            index += 1  # update index

    initial_concentration = np.array(initial_concentration,dtype=np.float64)    # convert concentration from list to array
    concentration = np.zeros((initial_concentration.shape[0],initial_concentration.shape[1],iters),dtype=np.float64)
    concentration[:,:,0] = initial_concentration.copy()
    sinking = np.array(sinking,dtype=np.float64)    # convert sinking from list to array

    return concentration, sinking, tracer_map, tracer_type


start = time.perf_counter()
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

# Import physical model
file = 'tests/npzd/physical_npzd.yaml'
# file = 'physical_bfm17_0d.yaml'
# file = 'physical_fasham_0d.yaml'
file_path = os.getcwd() + '/' + file
physical = import_physical_model(file_path)

file = 'tests/npzd/npzd.yaml'
# file = 'bfm17_0d.yaml'
# file = 'fasham_0d.yaml'
file_path = os.getcwd() + '/' + file
base_element, reactions, tracers = import_bgc_model(file_path, physical)

concentration, sinking, tracer_map, tracer_type = create_function_inputs(physical["simulation"]["iters"],tracers)


iters = physical["simulation"]["iters"]                             # iterations needed
dt = physical["simulation"]["dt"]                                   # time step [s]
num_layers = physical["water_column"]["num_layers"]                 # number of layers in water column [-]
column_depth = physical["water_column"]["column_depth"]             # water column depth [m]
z = physical["vertical_grid"]["z"]                                  # vertical grid [m]
dz = physical["vertical_grid"]["dz"]                                # vertical spacing [m]

configuration = physical["simulation"]["configuration"]
forcing = physical["environment"]["forcing"]
forcing_data = physical["environment"]["forcing_data"]
time_array = physical["simulation"]["time"]

# Initialize arrays
temperature = np.zeros(num_layers,dtype=np.float64)
salinity = np.zeros(num_layers,dtype=np.float64)
mixed_layer_depth = np.zeros(num_layers,dtype=np.float64)
surfacer_PAR = 0.
wind = np.zeros(num_layers,dtype=np.float64)

# Counters
day = 0
month = 1
print('month = ', month)
timesteps_per_day = 86400./dt

for iter in range(0,iters-1):
    if (iter != 0) & ((iter+1) % timesteps_per_day == 0): # take average at the end of day
        day += 1
        if (day != 0) & ((day+1) % 30 == 0):
            month += 1
            print('month = ', month)

    # Clear previous rates
    d_dt = np.zeros_like(concentration[...,iter])

    # Calculate physical variables at current time
    if forcing == "constant":
        temperature[0] = forcing_data["temperature"]
        salinity[0] = forcing_data["salinity"]
        mixed_layer_depth[0] = forcing_data["mld"]
        surfacer_PAR = forcing_data["sunlight"]
        wind[0] = forcing_data["wind"]
    elif forcing == "seasonal":
        temperature[0] = get_temperature(time_array[iter], forcing_data["winter_temp"], forcing_data["summer_temp"], forcing_data["temp_excursion"])
        salinity[0] = get_salinity(time_array[iter], forcing_data["winter_salt"], forcing_data["summer_salt"])
        mixed_layer_depth[0] = get_mixed_layer_depth(time_array[iter],forcing_data["winter_mld"], forcing_data["summer_mld"])
        surfacer_PAR = get_sunlight(time_array[iter],forcing_data["winter_sun"], forcing_data["summer_sun"], physical["environment"]["latitude"])
        wind[0] = get_wind(time_array[iter], forcing_data["winter_wind"], forcing_data["summer_wind"])
    density = density_profile(configuration, num_layers, column_depth/2, 0., temperature, salinity)     # Calculate density in center of cell (column_depth/2)

    # Calculate bgc rates
    d_dt = bgc_rate_eqns(iter, configuration, base_element, concentration[...,iter], d_dt, physical["environment"]["light_attenuation_water"], temperature, salinity, density, z, dz, surfacer_PAR, wind, tracer_map, tracer_type, tracers, sinking)

    # Update concentrations, set minimum of zero
    concentration[...,iter+1] = np.maximum(np.zeros_like(concentration[...,iter]), concentration[...,iter] + (dt * d_dt))

# ----------------------------------------------------------------------------------------------------
# Write outputs to npz file
# ----------------------------------------------------------------------------------------------------
npp_exists = False  # initialize writing of npp
for trac in tracers:
    if tracers[trac].type == "phytoplankton":
        if not npp_exists:  # first phytoplankton group
            npp = tracers[trac].npp
            npp_exists = True   # npp now exists, update to True to append with npp from later phytoplankton groups
        else:   # subsequent phytoplankton groups
            npp += tracers[trac].npp

conc_daily, conc_monthly = average(concentration,physical,'concentration')
np.savez('concentration_npzd_0908.npz',daily=conc_daily,monthly=conc_monthly)
if npp_exists:
    npp_daily, npp_monthly = average(npp,physical,'npp')
    np.savez('npp_npzd_0908.npz',daily=npp_daily,monthly=npp_monthly)

np.savez('tracer_indices_npzd_0908.npz',**tracer_map)

# ----------------------------------------------------------------------------------------------------
# Simulation complete
# ----------------------------------------------------------------------------------------------------
print('Simulation complete.')
elapsed = time.perf_counter() - start
hours, remainder = divmod(elapsed, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}")
