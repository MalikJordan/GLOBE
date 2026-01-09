import numpy as np
import os
import sys
import yaml
from setup.bacteria import Bacteria
from setup.detritus import Detritus
from setup.inorganic import Inorganic
from setup.phytoplankton import Phytoplankton
from setup.zooplankton import Zooplankton

# from setup.setup_solveivp import Bacteria, Detritus, Inorganic, Phytoplankton, Zooplankton


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
    iters = physical["simulation"]["days"] * 86400 / physical["simulation"]["timestep"]     
    iters = int(np.ceil(iters))

    # Create time array (needed for calculating incident angle for PAR)
    physical["simulation"]["time"] = np.linspace(0,iters * physical["simulation"]["timestep"], iters + 1)   # +1 for zero-indexing
    physical["simulation"]["iters"] = iters + 1 
    
    return physical



# def import_model(file_path):

#     with open(file_path, 'r') as f:
#         model_info = yaml.full_load(f)
#         base_element = model_info["base_element"]
#         model = model_info["tracers"]
#         parameters = model_info["parameters"]
#         reactions = model_info["reactions"]

#     # ----------------------------------------------------------------------------------------------------
#     # Check base element
#     # ----------------------------------------------------------------------------------------------------
#     if base_element not in ['c','n','p']:
#         sys.exit("'" + base_element + "' not accepted as base element. Check documentation and edit input file.")

#     # ----------------------------------------------------------------------------------------------------
#     # Update coordinate system
#     # ----------------------------------------------------------------------------------------------------
#     parameters["water_column"] = coordinate_system(parameters["water_column"])

#     # ----------------------------------------------------------------------------------------------------
#     # Setup time array
#     # ----------------------------------------------------------------------------------------------------
#     # Calculate number of iterations needed in simulation
#     iters = parameters["simulation"]["num_days"] * 86400 / parameters["simulation"]["timestep"]     
#     iters = int(np.ceil(iters))

#     # Create time array (needed for calculating incident angle for PAR)
#     parameters["simulation"]["time"] = np.linspace(0,iters * parameters["simulation"]["timestep"], iters + 1)   # +1 for zero-indexing
#     parameters["simulation"]["iters"] = iters + 1 
#     # ----------------------------------------------------------------------------------------------------
#     # Read model tracers
#     # ----------------------------------------------------------------------------------------------------
#     tracers = {}
#     for key in model:
#         if model[key]["type"] == "bacteria":
#             tracers[key] = Bacteria(key, parameters["simulation"]["iters"], reactions, **model[key])
#         elif model[key]["type"] == "detritus":
#             tracers[key] = Detritus(key, parameters["simulation"]["iters"], reactions, **model[key])
#         elif model[key]["type"] == "inorganic":
#             tracers[key] = Inorganic(key, parameters["simulation"]["iters"], reactions, **model[key])
#         elif model[key]["type"] == "phytoplankton":
#             tracers[key] = Phytoplankton(key, parameters["simulation"]["iters"], reactions, **model[key])
#         elif model[key]["type"] == "zooplankton":
#             tracers[key] = Zooplankton(key, parameters["simulation"]["iters"], reactions, **model[key])
#         else:
#             sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
            
#     # ----------------------------------------------------------------------------------------------------
#     # Add necessary components to Phytoplankton and Zooplankton groups
#     # ----------------------------------------------------------------------------------------------------
#     for key in reactions:
#         if key["type"] == "grazing":    # Add prey to zooplankton (used in rate calculations to determine sum of grazing rates)
#             tracers[list(key["produced"].keys())[0]].add_prey(list(key["consumed"].keys())[0])
#         if key["type"] == "uptake":     # Add nutrient to bacteria and phytoplankton (used in rate calculations for nutrient limitation)
#             # Add half saturation constant if used in calculations
#             tracers[list(key["produced"].keys())[0]].add_nutrient(list(key["consumed"].keys())[0])

#     return base_element, parameters, reactions, tracers


def initialize_pom(inputs):
    """
    Description: Opens forcing files reading the paths specified in the pom_input namelist.

    :return: data arrays for wind stress, surface salinity, solar radiation, inorganic
             suspended matter, salinity and temperature vertical profiles, general circulation
             for w velocity, intermediate eddy velocities, salinity and temperature initial
             conditions, heat flux loss, and surface and bottom nutrients
    """

    inputs["input_files"]
    # vertical_layers = inputs["water_column"]["num_layers"]
    # Length of input arrays
    array_length = 13   # months (D-J-F-M-A-M-J-J-A-S-O-N-D)

    # Wind speed
    wind_speed_data = np.fromfile(inputs["input_files"]["wind_stress"])
    wind_speed_zonal   = np.zeros(array_length)
    wind_speed_meridional   = np.zeros(array_length)
    for i in range(0,array_length):
        wind_speed_zonal[i] = wind_speed_data[2*i + 0]
        wind_speed_meridional[i] = wind_speed_data[2*i + 1]

    # SOlar radiation, shortwave radiation, and heat flux
    solar_radiation = np.fromfile(inputs["input_files"]["shortwave_solar_radiation"])
    heat_flux_loss_data = np.fromfile(inputs["input_files"]["heat_flux_loss"])
    shortwave_radiation = np.zeros(array_length)
    surface_heat_flux = np.zeros(array_length)
    kinetic_energy_loss = np.zeros(array_length)
    for i in range(0,array_length):
        shortwave_radiation[i]  = heat_flux_loss_data[3*i + 0]
        surface_heat_flux[i] = heat_flux_loss_data[3*i + 1]
        kinetic_energy_loss[i]  = heat_flux_loss_data[3*i + 2]

    # Inorganic suspended matter
    inorganic_suspended_matter_data = np.fromfile(inputs["input_files"]["shortwave_solar_radiation"])
    inorganic_suspended_matter   = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            inorganic_suspended_matter[x,i] = inorganic_suspended_matter_data[inputs["water_column"]["num_layers"] * i + x]

    # Surface salinity
    surface_salinity = np.fromfile(inputs["input_files"]["surface_salinity"])

    # Salinity climatology (diagnostic mode)
    salinity_vertical_profile_data = np.fromfile(inputs["input_files"]["salinity"])
    salinity_climatology = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            salinity_climatology[x,i] = salinity_vertical_profile_data[inputs["water_column"]["num_layers"] * i + x]

    # Salinity IC
    salinity = np.fromfile(inputs["input_files"]["salinity_IC"])

    # Temperature climatology (diagnostic mode)
    temperature_vertical_profile_data = np.fromfile(inputs["input_files"]["temperature"])
    temperature_climatology = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            temperature_climatology[x,i] = temperature_vertical_profile_data[inputs["water_column"]["num_layers"] * i + x]

    # Temperature IC
    temperature = np.fromfile(inputs["input_files"]["temperature_IC"])

    # General circulation w velocity climatology
    general_circulation_w_velocity_data = np.fromfile(inputs["input_files"]["w_velocity"])
    w_velocity_climatology  = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            w_velocity_climatology[x,i] = general_circulation_w_velocity_data[inputs["water_column"]["num_layers"] * i + x]

    # Intermittant eddy w velocity 1
    intermediate_eddy_w_velocity_1_data = np.fromfile(inputs["input_files"]["eddy_w_velocity_1"])
    w_eddy_velocity_1  = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            w_eddy_velocity_1[x,i] = intermediate_eddy_w_velocity_1_data[inputs["water_column"]["num_layers"] * i + x]


    # Intermittant eddy w velocity 2
    intermediate_eddy_w_velocity_2_data = np.fromfile(inputs["input_files"]["eddy_w_velocity_2"])
    w_eddy_velocity_2  = np.zeros((inputs["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, inputs["water_column"]["num_layers"]):
            w_eddy_velocity_1[x,i] = intermediate_eddy_w_velocity_1_data[inputs["water_column"]["num_layers"] * i + x]

    # Surface nutrients
    surface_nutrients_data  = np.fromfile(inputs["input_files"]["surface_nutrients"])
    NO3_s1  = np.zeros(array_length)
    NH4_s1  = np.zeros(array_length)
    PO4_s1  = np.zeros(array_length)
    SIO4_s1 = np.zeros(array_length)
    for i in range(0,array_length):
        NO3_s1[i]  = surface_nutrients_data[4*i + 0]
        NH4_s1[i]  = surface_nutrients_data[4*i + 1]
        PO4_s1[i]  = surface_nutrients_data[4*i + 2]
        SIO4_s1[i] = surface_nutrients_data[4*i + 3]

    # Bottom nutrients
    bottom_nutrients_data = np.fromfile(inputs["input_files"]["bottom_nutrients"])
    O2_b1   = np.zeros(array_length)
    NO3_b1  = np.zeros(array_length)
    PO4_b1  = np.zeros(array_length)
    PON_b1  = np.zeros(array_length)
    for i in range(0,array_length):
        O2_b1[i]  = bottom_nutrients_data[4*i + 0]
        NO3_b1[i] = bottom_nutrients_data[4*i + 1]
        PO4_b1[i] = bottom_nutrients_data[4*i + 2]
        PON_b1[i] = bottom_nutrients_data[4*i + 3]

    

def temperature_salinity_initial_conditions(inputs):

    # Salinity
    salinity = np.fromfile(inputs["input_files"]["salinity_IC"])

    # Temperature
    temperature = np.fromfile(inputs["input_files"]["temperature_IC"])




