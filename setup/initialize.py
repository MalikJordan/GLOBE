import numpy as np
import sys
import yaml
from setup.bacteria import Bacteria
from setup.detritus import Detritus
from setup.inorganic import Inorganic
from setup.phytoplankton import Phytoplankton
from setup.zooplankton import Zooplankton


def coordinate_system(parameters):

    parameters["dz"] = parameters["column_depth"] / parameters["num_layers"]
    parameters["z"] = np.linspace(parameters["dz"]/2, parameters["column_depth"] - parameters["dz"]/2, parameters["num_layers"])

    return parameters


def import_model(file_path):

    with open(file_path, 'r') as f:
        model_info = yaml.full_load(f)
        base_element = model_info["base_element"]
        model = model_info["tracers"]
        parameters = model_info["parameters"]
        reactions = model_info["reactions"]

    # ----------------------------------------------------------------------------------------------------
    # Check base element
    # ----------------------------------------------------------------------------------------------------
    if base_element not in ['c','n','p']:
        sys.exit("'" + base_element + "' not accepted as base element. Check documentation and edit input file.")

    # ----------------------------------------------------------------------------------------------------
    # Update coordinate system
    # ----------------------------------------------------------------------------------------------------
    parameters["water_column"] = coordinate_system(parameters["water_column"])

    # ----------------------------------------------------------------------------------------------------
    # Setup time array
    # ----------------------------------------------------------------------------------------------------
    # Calculate number of iterations needed in simulation
    parameters["simulation"]["iters"] = parameters["simulation"]["num_days"] * 86400 / parameters["simulation"]["timestep"]     
    parameters["simulation"]["iters"] = int(np.ceil(parameters["simulation"]["iters"]))

    # Create time array (needed for calculating incident angle for PAR)
    parameters["simulation"]["time"] = np.linspace(0,parameters["simulation"]["iters"] * parameters["simulation"]["timestep"], parameters["simulation"]["iters"] + 1)   # +1 for zero-indexing
    # ----------------------------------------------------------------------------------------------------
    # Read model tracers
    # ----------------------------------------------------------------------------------------------------
    tracers = {}
    for key in model:
        if model[key]["type"] == "bacteria":
            # tracers[key] = Bacteria(key, model[key]["composition"], parameters["simulation"]["iters"], model[key]["long_name"], model[key]["parameters"], reactions, model[key]["type"])
            tracers[key] = Bacteria(parameters["simulation"]["iters"], reactions, **model[key])
        elif model[key]["type"] == "detritus":
            # tracers[key] = Detritus(key, model[key]["composition"], parameters["simulation"]["iters"], model[key]["long_name"], reactions, model[key]["type"])
            tracers[key] = Detritus(parameters["simulation"]["iters"], reactions, **model[key])
        elif model[key]["type"] == "inorganic":
            # tracers[key] = Inorganic(key, model[key]["composition"], parameters["simulation"]["iters"], model[key]["long_name"], reactions, model[key]["type"])
            tracers[key] = Inorganic(parameters["simulation"]["iters"], reactions, **model[key])
        elif model[key]["type"] == "phytoplankton":
            # tracers[key] = Phytoplankton(key, model[key]["composition"], parameters["simulation"]["iters"], model[key]["long_name"], model[key]["parameters"], reactions, model[key]["type"])
            tracers[key] = Phytoplankton(parameters["simulation"]["iters"], reactions, **model[key])
        elif model[key]["type"] == "zooplankton":
            # tracers[key] = Zooplankton(key, model[key]["composition"], parameters["simulation"]["iters"], model[key]["long_name"], model[key]["parameters"], reactions, model[key]["type"])
            tracers[key] = Zooplankton(parameters["simulation"]["iters"], reactions, **model[key])
        else:
            sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
            
    # ----------------------------------------------------------------------------------------------------
    # Add necessary components to Phytoplankton and Zooplankton groups
    # ----------------------------------------------------------------------------------------------------
    for key in reactions:
        if key["type"] == "grazing":    # Add prey to zooplankton (used in rate calculations to determine sum of grazing rates)
            # produced = list(key["produced"].keys())[0]
            # consumed = list(key["consumed"].keys())[0]
            tracers[list(key["produced"].keys())[0]].add_prey(list(key["consumed"].keys())[0])
        if key["type"] == "uptake":     # Add nutrient to bacteria and phytoplankton (used in rate calculations for nutrient limitation)
            # produced = list(key["produced"].keys())[0]
            # consumed = list(key["consumed"].keys())[0]
            tracers[list(key["produced"].keys())[0]].add_nutrient(list(key["consumed"].keys())[0],key["parameters"]["half_sat_nutrient"])

    return base_element, parameters, reactions, tracers
