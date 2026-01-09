import numpy as np
import os
import sys
import yaml
from setup.bacteria import Bacteria
from setup.detritus import Detritus
from setup.inorganic import Inorganic
from setup.phytoplankton import Phytoplankton
from setup.zooplankton import Zooplankton


def coordinate_system(configuration, parameters):

    # Create empty dictionary for vertical grid parameters
    vertical_grid = {}

    # Create coordinate system
    if configuration == "0d":
        vertical_grid["dz"] = parameters["column_depth"] / parameters["num_layers"]
        vertical_grid["z"] = np.linspace(parameters["dz"]/2, parameters["column_depth"] - parameters["dz"]/2, parameters["num_layers"])

    elif configuration == "1d":
        # Initialize vertical coordinate system arrays
        l = np.ones(parameters["num_layers"])       # length scale
        z = np.zeros(parameters["num_layers"])      # vertical coordinates
        zz = np.zeros(parameters["num_layers"])     # staggered vertical coordinates
        dz = np.zeros(parameters["num_layers"])     # vertical spacing
        dzz = np.zeros(parameters["num_layers"])    # staggered vertical spacing
        dzr = np.zeros(parameters["num_layers"])    # reciprocal of vertical spacing

        # Calculate initial spacing
        surface_logspace_layers = parameters["surf_log"] - 2.
        bottom_logspace_layers = parameters["num_layers"] - parameters["bot_log"] - 1.

        layers = (parameters["bot_log"] - parameters["surf_log"]) + 4.
        
        initial_spacing = 2. / layers / np.exp(0.693147 * (surface_logspace_layers))

        dzz[0] = -0.5 * initial_spacing

        # Set vertical coordinates
        for i in range(1, int(parameters["surf_log"])-1):
            z[i-1] = -initial_spacing * 2**(i-2)
            zz[i-1] = -initial_spacing * 2**(i-1.5)

        for i in range(int(parameters["surf_log"])-1, parameters["num_layers"]+1):
            z[i-1] = -(i - surface_logspace_layers) / layers
            zz[i-1] = -(i - surface_logspace_layers + 0.5) / layers
        
        # Set vertical spacing
        dz[:-1] = z[:-1] - z[1:]
        dzz[:-1] = zz[:-1] - zz[1:]

        dz[-1] = 1.E-06     # Small value to avoid division by zero for dzr
        dzr = 1. / dz       # Take reciprocal
        dz[-1] = 0.         # Correct dz value

        # Set length scale
        l[0] = 0.
        l[-1] = 0.

        vertical_grid["l"] = l
        vertical_grid["z"] = z
        vertical_grid["zz"] = zz
        vertical_grid["dz"] = dz
        vertical_grid["dzz"] = dzz
        vertical_grid["dzr"] = dzr

    return vertical_grid


def import_bgc_model(file_path, physical):

    # Open file containing bgc data
    with open(file_path, 'r') as f:
        model_info = yaml.full_load(f)
        base_element = model_info["base_element"]
        model = model_info["tracers"]
        reactions = model_info["reactions"]

    # ----------------------------------------------------------------------------------------------------
    # Read model tracers
    # ----------------------------------------------------------------------------------------------------
    tracers = {}
    for key in model:
        if model[key]["type"] == "bacteria":
            tracers[key] = Bacteria(key, base_element, physical["simulation"]["iters"], physical["water_column"]["num_layers"], reactions, **model[key])
        elif model[key]["type"] == "detritus":
            tracers[key] = Detritus(key, base_element, physical["simulation"]["iters"], physical["water_column"]["num_layers"], reactions, **model[key])
        elif model[key]["type"] == "inorganic":
            tracers[key] = Inorganic(key, physical["simulation"]["iters"], physical["water_column"]["num_layers"], reactions, **model[key])
        elif model[key]["type"] == "phytoplankton":
            tracers[key] = Phytoplankton(key, base_element, physical["simulation"]["iters"], physical["water_column"]["num_layers"], reactions, **model[key])
        elif model[key]["type"] == "zooplankton":
            tracers[key] = Zooplankton(key, base_element, physical["simulation"]["iters"], physical["water_column"]["num_layers"], reactions, **model[key])
        else:
            sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
    
    # tracers = {}
    # for key in model:
    #     if model[key]["type"] == "bacteria":
    #         tracers[key] = Bacteria(key, physical["simulation"]["iters"], reactions, **model[key])
    #     elif model[key]["type"] == "detritus":
    #         tracers[key] = Detritus(key, physical["simulation"]["iters"], reactions, **model[key])
    #     elif model[key]["type"] == "inorganic":
    #         tracers[key] = Inorganic(key, physical["simulation"]["iters"], reactions, **model[key])
    #     elif model[key]["type"] == "phytoplankton":
    #         tracers[key] = Phytoplankton(key, physical["simulation"]["iters"], reactions, **model[key])
    #     elif model[key]["type"] == "zooplankton":
    #         tracers[key] = Zooplankton(key, physical["simulation"]["iters"], reactions, **model[key])
    #     else:
    #         sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
            
    # ----------------------------------------------------------------------------------------------------
    # Add necessary components to Phytoplankton and Zooplankton groups
    # ----------------------------------------------------------------------------------------------------
    for key in reactions:
        if key["type"] == "grazing":    # Add prey to zooplankton (used in rate calculations to determine sum of grazing rates)
            tracers[list(key["produced"].keys())[0]].add_prey(list(key["consumed"].keys())[0])
        if key["type"] == "uptake":     # Add nutrient to bacteria and phytoplankton (used in rate calculations for nutrient limitation)
            # Add half saturation constant if used in calculations
            tracers[list(key["produced"].keys())[0]].add_nutrient(list(key["consumed"].keys())[0])

    return base_element, reactions, tracers


def import_physical_model(file_path):
    
    # Open file containing physical data
    with open(file_path, 'r') as f:     physical = yaml.full_load(f)
    
    # ----------------------------------------------------------------------------------------------------
    # Initialize coordinate system
    # ----------------------------------------------------------------------------------------------------
    physical["vertical_grid"] = coordinate_system(physical["simulation"]["configuration"],physical["water_column"])

    # ----------------------------------------------------------------------------------------------------
    # Setup time array
    # ----------------------------------------------------------------------------------------------------
    # Calculate number of iterations needed in simulation
    physical["simulation"]["sec_per_day"] = 86400
    iters = physical["simulation"]["days"] * physical["simulation"]["sec_per_day"] / physical["simulation"]["dt"]     
    iters = int(np.ceil(iters))

    # Create time array (needed for calculating incident angle for PAR)
    physical["simulation"]["time"] = np.linspace(0,iters * physical["simulation"]["dt"], iters + 1)   # +1 for zero-indexing
    physical["simulation"]["iters"] = iters + 1 

    # Other variables for future calculations
    physical["simulation"]["dt2"] = 2. * physical["simulation"]["dt"]   # twice the timestep
    physical["simulation"]["day_per_sec"] = 1. / 86400
    physical["simulation"]["months"] = int(np.ceil(physical["simulation"]["days"]/30))
    
    return physical
