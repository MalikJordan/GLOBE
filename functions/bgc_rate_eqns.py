import numpy as np
import os
import yaml
from functions import rates, seasonal_cycling
from functions.other_functions import concentration_ratio
from setup.initialize import coordinate_system


def bgc_rate_eqns(iter, base_element, parameters, tracers):
# def bgc_rate_eqns(iter):
    
#     from globe import base_element, parameters, tracers
    # Extract model paramters
    environmental_parameters = parameters["environment"]
    simulation_parameters = parameters["simulation"]
    water_column_parameters = parameters["water_column"]
    water_column_parameters = coordinate_system(parameters["water_column"])

    # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
    temperature = seasonal_cycling.get_temperature(simulation_parameters["time"][iter], environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    surface_PAR = seasonal_cycling.get_sunlight(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"],environmental_parameters["latitude"])
    # surface_PAR = seasonal_cycling.get_irrad(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(simulation_parameters["time"][iter], environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])
    salinity = seasonal_cycling.get_salinity(simulation_parameters["time"][iter], environmental_parameters["winter_salt"], environmental_parameters["summer_salt"])
    wind = seasonal_cycling.get_wind(simulation_parameters["time"][iter], environmental_parameters["winter_wind"], environmental_parameters["summer_wind"])

    # Clear previous rates
    for key in tracers:
        tracers[key].d_dt = np.zeros_like(tracers[key].d_dt)
    
    if iter == 360:
        x=1
    
    # Calculate bgc rates
    for key in tracers:
        if tracers[key].type == "bacteria":
            pass
        elif tracers[key].type == "detritus":
            tracers[key].detritus(iter, base_element, tracers)
        elif tracers[key].type == "inorganic":
            tracers[key].inorg(iter, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, salinity, wind, tracers)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(iter, base_element, environmental_parameters["base_temp"], environmental_parameters["light_attenuation_water"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, tracers)
            # fI, irrad, nut_lim = tracers[key].phyto(iter, base_element, environmental_parameters["base_temp"], environmental_parameters["light_attenuation_water"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, tracers)
            # fI, irrad = tracers[key].phyto(iter, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, tracers)
        elif tracers[key].type == "zooplankton":
            tracers[key].zoo(iter, base_element, environmental_parameters["base_temp"], temperature, tracers)
            # rsp = tracers[key].zoo(iter, base_element, environmental_parameters["base_temp"], temperature, tracers)
        
    # Apply rates to tracer concentrations
    for key in tracers:
        tracers[key].conc[...,iter+1] = tracers[key].conc[...,iter] + simulation_parameters["timestep"] * tracers[key].d_dt / 86400  # Convert rates from 1/d to 1/s
    
    # Set minimum of zero for tracer concentrations
    for key in tracers:
        tracers[key].conc[...,iter+1] = np.maximum(1E-20*np.ones_like(tracers[key].conc[...,iter+1]),tracers[key].conc[...,iter+1])
    
    # Update concentration ratios
    for key in tracers:
        if tracers[key].type in ["bacteria", "detritus","phytoplankton","zooplankton"]:
            # Get index of base element
            index = tracers[key].composition.index(base_element)

            # Calculate concentration ratios
            concentration_ratio(iter+1, index, tracers[key])
    # return fI, irrad, nut_lim


def bgc_rate_eqns_solveivp(time, concentration):
    from test_globe import base_element, parameters, tracers, tracer_indices
    # Extract model paramters
    environmental_parameters = parameters["environment"]
    simulation_parameters = parameters["simulation"]
    water_column_parameters = parameters["water_column"]
    water_column_parameters = coordinate_system(parameters["water_column"])

    # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
    temperature = seasonal_cycling.get_temperature(time, environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    # surface_PAR = seasonal_cycling.get_sunlight(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    surface_PAR = seasonal_cycling.get_irrad(time, environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(time, environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])
    salinity = seasonal_cycling.get_salinity(time, environmental_parameters["winter_salt"], environmental_parameters["summer_salt"])
    wind = seasonal_cycling.get_wind(time, environmental_parameters["winter_wind"], environmental_parameters["summer_wind"])

    # temperature = seasonal_cycling.get_temperature(simulation_parameters["time"][iter], environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    # # surface_PAR = seasonal_cycling.get_sunlight(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    # surface_PAR = seasonal_cycling.get_irrad(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    # mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(simulation_parameters["time"][iter], environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])
    # salinity = seasonal_cycling.get_salinity(simulation_parameters["time"][iter], environmental_parameters["winter_salt"], environmental_parameters["summer_salt"])
    # wind = seasonal_cycling.get_wind(simulation_parameters["time"][iter], environmental_parameters["winter_wind"], environmental_parameters["summer_wind"])
    
    # Clip concentration matrix
    for i in range(0,len(concentration)):
        concentration[i] = np.maximum(1E-20,concentration[i])

    # Clear previous rates
    for key in tracers:
        tracers[key].d_dt = np.zeros_like(tracers[key].d_dt)
    
    # Calculate bgc rates
    for key in tracers:
        if tracers[key].type == "bacteria":
            pass
        elif tracers[key].type == "detritus":
            tracers[key].detritus(concentration, tracer_indices, base_element, tracers)
        elif tracers[key].type == "inorganic":
            tracers[key].inorg(concentration, tracer_indices, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, salinity, wind, tracers)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(concentration, tracer_indices, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, tracers)
            # fI, irrad = tracers[key].phyto(iter, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, tracers)
        elif tracers[key].type == "zooplankton":
            tracers[key].zoo(concentration, tracer_indices, base_element, environmental_parameters["base_temp"], temperature, tracers)
            # rsp = tracers[key].zoo(iter, base_element, environmental_parameters["base_temp"], temperature, tracers)

    # Update rate of change matrix
    dt = np.zeros_like(concentration)
    for t in tracers:
        indices = tracer_indices[t]
        for i in range(len(indices)):
            dt[indices[i]] = tracers[t].d_dt[i,...]

    return dt