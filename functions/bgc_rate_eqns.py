import json
import numpy as np
import os
from functions import rates, seasonal_cycling
from setup.initialize import coordinate_system

# Open json file containing reactions and model parameters
file_path = os.getcwd() + "/npzd-riley.json"
# file_path = os.getcwd() + "/model_description.json"
with open(file_path) as read_model:
    model_info = json.load(read_model)
    model = model_info["tracers"]
    parameters = model_info["parameters"]
    reactions = model_info["reactions"]

# Extract model parameters
environmental_parameters = parameters["environment"]
simulation_parameters = parameters["simulation"]
water_column_parameters = parameters["water_column"]
water_column_parameters = coordinate_system(water_column_parameters)


def bgc_rate_eqns(time, conc, tracers):

    bgc_rates = np.zeros_like(conc)

    # Seasonal cycling for temperature, salinity, radiation, and mixed layer depth
    temperature = seasonal_cycling.get_temperature(time, environmental_parameters["winter_temp"], environmental_parameters["summer_temp"])
    surface_PAR = seasonal_cycling.get_sunlight(time, environmental_parameters["winter_sun"], environmental_parameters["summer_sun"])
    mixed_layer_depth = seasonal_cycling.get_mixed_layer_depth(time, environmental_parameters["winter_mld"], environmental_parameters["summer_mld"])

    # If "Grazing" is listed as a reaction, create dictionary to store a list of grazing rates for all zooplankton groups
    grazing_rates = None 

    for reac in reactions:
        if reac["type"] == "grazing":
            if grazing_rates is None:
                grazing_rates = {}
    if grazing_rates is not None:
        for key in model["zooplankton"]:
            grazing_rates[key] = []

    for reac in reactions:
        if reac["type"] == "grazing":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate grazing rate
            d_dt_consumed, d_dt_produced = rates.grazing(reac["parameters"], consumed, produced)

            # Append grazing list
            grazing_rates[reac["produced"]].append(d_dt_consumed)

            # Update rates for consumed (-) and produced (+) tracers
            bgc_rates[tracers[reac["produced"]].index] += d_dt_produced
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt_consumed

        if reac["type"] == "egestion":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate egestion rate
            d_dt = rates.egestion(reac["parameters"], grazing_rates[reac["consumed"]])

            # Update rates for consumed (-) and produced (+) tracers
            bgc_rates[tracers[reac["produced"]].index] += d_dt
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt

        if reac["type"] == "excretion":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate egestion rate
            d_dt = rates.excretion(reac["parameters"], grazing_rates[reac["consumed"]])

            # Update rates for consumed (-) and produced (+) tracers
            bgc_rates[tracers[reac["produced"]].index] += d_dt
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt

        if reac["type"] == "mortality":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate egestion rate
            d_dt = rates.mortality(reac["parameters"], consumed)

            bgc_rates[tracers[reac["produced"]].index] += d_dt
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt
            # Update rates for consumed (-) and produced (+) tracers
            # if "zooplankton" in model and reac["consumed"] in model["zooplankton"]: # Only linear term (natural) for zooplankton is allocated to detrital pool
            #     bgc_rates[tracers[reac["produced"]].index] += d_dt_natural  # Add natural mortality to detrital pool
            #     bgc_rates[tracers[reac["consumed"]].index] -= d_dt_total    # Subtract total mortality from zooplankton concentration
            # else:
            #     bgc_rates[tracers[reac["produced"]].index] += d_dt_total
            #     bgc_rates[tracers[reac["consumed"]].index] -= d_dt_total

        if reac["type"] == "remineralization":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate remineralization rate
            d_dt = rates.remineralization(reac["parameters"], consumed)

            # Update rates for consumed (-) and produced (+) tracers
            bgc_rates[tracers[reac["produced"]].index] += d_dt
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt

        if reac["type"] == "uptake":
            # Extract tracer concentrations
            consumed = conc[tracers[reac["consumed"]].index]
            produced = conc[tracers[reac["produced"]].index]

            # Calculate growth rate
            d_dt = rates.uptake(reac["parameters"], environmental_parameters["base_temp"], water_column_parameters["z"], mixed_layer_depth, consumed, produced, surface_PAR, temperature)

            # Update rates for consumed (-) and produced (+) tracers
            bgc_rates[tracers[reac["produced"]].index] += d_dt
            bgc_rates[tracers[reac["consumed"]].index] -= d_dt

    bgc_rates /= 86400  # Convert rates from [1/d] to [1/s]

    return bgc_rates
