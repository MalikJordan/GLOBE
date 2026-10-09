import copy
import os
import yaml
from optimization.build import extract_parameters


def optimize_bgc_model(base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical):

    if physical["simulation"]["reduce"] == True:    model_file_path = os.getcwd() + '/reduced_model.yaml' 
    else:   model_file_path = os.getcwd() + '/model.yaml'

    # Read bgc input file
    with open(model_file_path, 'r') as f:
        model_info = yaml.full_load(f)
        model = model_info["tracers"]

    # Begin optimization
    if physical["simulation"]["optimize"] == True:
        # Read dakota input file
        dakota_file_path = os.getcwd() + '/dakota.yaml'
        with open(dakota_file_path, 'r') as f:
            dakota_configuration = yaml.full_load(f)

        # Extract parameters and calculate bounds
        parameters, numeric_parameters, bounds = extract_parameters(model, dakota_configuration)

        # Generate and run Dakota Morris screening on numeric parameter set
        

    return base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical