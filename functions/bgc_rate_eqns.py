import numpy as np
import os
import yaml
from functions import rates, seasonal_cycling
from setup.initialize import coordinate_system


def bgc_rate_eqns(time, file_path, tracers):
    
    # ----------------------------------------------------------------------------------------------------
    # Extract model paramters
    # ----------------------------------------------------------------------------------------------------
    with open(file_path, 'r') as f:
        model_info = yaml.full_load(f)
        base_element = model_info["base_element"]
        model = model_info["tracers"]
        parameters = model_info["parameters"]
        reactions = model_info["reactions"]

    environmental_parameters = parameters["environment"]
    simulation_parameters = parameters["simulation"]
    water_column_parameters = parameters["water_column"]
    water_column_parameters = coordinate_system(parameters["water_column"])

    # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
    temperature = seasonal_cycling.get_temperature(time, environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    surface_PAR = seasonal_cycling.get_sunlight(time, environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(time, environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])

    # Clear previous rates
    for key in tracers:
        tracers[key].d_dt = np.zeros_like(tracers[key].d_dt)
    
    # Calculate bgc rates
    for key in tracers:
        if tracers[key].type == "bacteria":
            pass
        elif tracers[key].type == "detritus":
            pass
            # tracers[key].detritus(base_element, reactions, tracers)
        elif tracers[key].type == "inorganic":
            pass
        elif tracers[key].type == "phytoplankton":
            pass
            # tracers[key].phyto(base_element, environmental_parameters["base_temp"], water_column_parameters["z"], mixed_layer_depth, reactions, surface_PAR, temperature, tracers)
        elif tracers[key].type == "zooplankton":
            tracers[key].zoo(base_element, environmental_parameters["base_temp"], temperature, tracers)
        
    # Apply rates to tracer concentrations
    for key in tracers:
        tracers[key].conc += simulation_parameters["timestep"] * tracers[key].d_dt / 86400  # Convert rates from 1/d to 1/s

    return 


# def bgc_rate_eqns(time, file_path, conc, tracers):

#     # ----------------------------------------------------------------------------------------------------
#     # Extract model paramters
#     # ----------------------------------------------------------------------------------------------------
#     with open(file_path, 'r') as f:
#         model_info = yaml.full_load(f)
#         base_element = model_info["base_element"]
#         model = model_info["tracers"]
#         parameters = model_info["parameters"]
#         reactions = model_info["reactions"]

#     environmental_parameters = parameters["environment"]
#     simulation_parameters = parameters["simulation"]
#     water_column_parameters = parameters["water_column"]
#     water_column_parameters = coordinate_system(parameters["water_column"])

#     bgc_rates = np.zeros_like(conc)

#     # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
#     temperature = seasonal_cycling.get_temperature(time, environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
#     surface_PAR = seasonal_cycling.get_sunlight(time, environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
#     mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(time, environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])

#     # If "Grazing" is listed as a reaction, create dictionary to store a list of grazing rates for all zooplankton groups
#     grazing_rates = None 

#     for reac in reactions:
#         if reac["type"] == "grazing":
#             if grazing_rates is None:
#                 grazing_rates = {}
#             break
#     if grazing_rates is not None:
#         for key in model["zooplankton"]:
#             grazing_rates[key] = []

#     for reac in reactions:
#         if reac["type"] == "grazing":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate grazing rate
#             d_dt_consumed, d_dt_produced = rates.grazing(reac["parameters"], consumed, produced)

#             # Append grazing list
#             grazing_rates[reac["produced"]].append(d_dt_consumed)

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt_produced
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt_consumed

#         if reac["type"] == "egestion":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate egestion rate
#             d_dt = rates.egestion(reac["parameters"], grazing_rates[reac["consumed"]])

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt

#         if reac["type"] == "excretion":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate egestion rate
#             d_dt = rates.excretion(reac["parameters"], grazing_rates[reac["consumed"]], consumed)

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt

#         if reac["type"] == "mortality":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate mortality rate
#             d_dt = rates.mortality(reac["parameters"], consumed)

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt

#         if reac["type"] == "remineralization":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate remineralization rate
#             d_dt = rates.remineralization(reac["parameters"], consumed)

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt

#         if reac["type"] == "uptake":
#             # Extract tracer concentrations
#             consumed = conc[tracers[reac["consumed"]].index]
#             produced = conc[tracers[reac["produced"]].index]

#             # Calculate growth rate
#             d_dt = rates.uptake(reac["parameters"], environmental_parameters["base_temp"], water_column_parameters["z"], mixed_layer_depth, consumed, produced, surface_PAR, temperature)

#             # Update rates for consumed (-) and produced (+) tracers
#             bgc_rates[tracers[reac["produced"]].index] += d_dt
#             bgc_rates[tracers[reac["consumed"]].index] -= d_dt

#     bgc_rates /= 86400  # Convert rates from [1/d] to [1/s]

#     return bgc_rates
