import numpy as np
import os
import yaml
from functions import rates, seasonal_cycling
from setup.initialize import coordinate_system


def bgc_rate_eqns(iter, base_element, parameters, tracers):
    
    # Extract model paramters
    environmental_parameters = parameters["environment"]
    simulation_parameters = parameters["simulation"]
    water_column_parameters = parameters["water_column"]
    water_column_parameters = coordinate_system(parameters["water_column"])

    # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
    temperature = seasonal_cycling.get_temperature(simulation_parameters["time"][iter], environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    surface_PAR = seasonal_cycling.get_sunlight(simulation_parameters["time"][iter], environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(simulation_parameters["time"][iter], environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])
    salinity = seasonal_cycling.get_salinity(simulation_parameters["time"][iter], environmental_parameters["winter_salt"], environmental_parameters["summer_salt"])
    wind = seasonal_cycling.get_wind(simulation_parameters["time"][iter], environmental_parameters["winter_wind"], environmental_parameters["summer_wind"])

    # Clear previous rates
    for key in tracers:
        tracers[key].d_dt = np.zeros_like(tracers[key].d_dt)
    
    # Calculate bgc rates
    for key in tracers:
        if tracers[key].type == "bacteria":
            pass
        elif tracers[key].type == "detritus":
            tracers[key].detritus(iter, base_element, tracers)
        elif tracers[key].type == "inorganic":
            tracers[key].inorg(iter, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], water_column_parameters["dz"], mixed_layer_depth, surface_PAR, temperature, salinity, wind, tracers)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(iter, base_element, environmental_parameters["base_temp"], water_column_parameters["z"], mixed_layer_depth, surface_PAR, temperature, tracers)
        elif tracers[key].type == "zooplankton":
            tracers[key].zoo(iter, base_element, environmental_parameters["base_temp"], temperature, tracers)
        
    # Apply rates to tracer concentrations
    for key in tracers:
        tracers[key].conc[...,iter+1] = tracers[key].conc[...,iter] + simulation_parameters["timestep"] * tracers[key].d_dt / 86400  # Convert rates from 1/d to 1/s
    
    # Set minimum of zero for tracer concentrations
    for key in tracers:
        tracers[key].conc[...,iter+1] = np.maximum(np.ones_like(tracers[key].conc[...,iter+1]),tracers[key].conc[...,iter+1])
        
    return