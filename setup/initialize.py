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
    # Read model tracers
    # ----------------------------------------------------------------------------------------------------
    tracers = {}
    for key in model:
        if model[key]["type"] == "bacteria":
            tracers[key] = Bacteria(model[key]["long_name"], model[key]["composition"], model[key]["type"])
        elif model[key]["type"] == "detritus":
            tracers[key] = Detritus(model[key]["long_name"], model[key]["composition"], model[key]["type"])
        elif model[key]["type"] == "inorganic":
            tracers[key] = Inorganic(model[key]["long_name"], model[key]["composition"], model[key]["type"])
        elif model[key]["type"] == "phytoplankton":
            tracers[key] = Phytoplankton(model[key]["long_name"], model[key]["composition"], model[key]["type"], model[key]["parameters"]["nutrient_limitation"])
        elif model[key]["type"] == "zooplankton":
            tracers[key] = Zooplankton(key, model[key]["long_name"], model[key]["composition"], model[key]["parameters"], model[key]["type"], reactions)
        # elif model[key]["type"] == "microzooplankton":
        #     tracers[key] = Zooplankton(model[key]["long_name"], model[key]["composition"], model[key]["type"])
        else:
            sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
            
    # ----------------------------------------------------------------------------------------------------
    # Add necessary components to Phytoplankton and Zooplankton groups
    # ----------------------------------------------------------------------------------------------------
    for key in reactions:
        if key["type"] == "grazing":    # Add prey to zooplankton (used in rate calculations to determine sum of grazing rates)
            # tracers[key["produced"]].add_prey(key["consumed"])
            produced = list(key["produced"].keys())[0]
            consumed = list(key["consumed"].keys())[0]
            tracers[produced].add_prey(consumed)
        if key["type"] == "uptake":     # Add nutrient to phytoplankton (used in rate calculations for nutrient limitation)
            # tracers[key["produced"]].add_nutrient(key["consumed"],key["parameters"]["half_sat_nutrient"])
            produced = list(key["produced"].keys())[0]
            consumed = list(key["consumed"].keys())[0]
            tracers[produced].add_nutrient(consumed,key["parameters"]["half_sat_nutrient"])

    return base_element, parameters, reactions, tracers
