import os
import yaml
from pathlib import Path
from setup.other_functions import create_function_inputs
from setup.initialize import import_bgc_model, import_physical_model

def initialize_model():
    """
    Definition:: handles the initialization of bgc model and physical/environmental parameters
    """

    print('Initializing model...')
    # physical_file_path = os.getcwd() + '/physical.yaml'
    physical_file_path = Path.cwd() / 'physical.yaml'
    physical = import_physical_model(physical_file_path)

    # model_file_path = os.getcwd() + '/model.yaml'
    model_file_path = Path.cwd() / 'model.yaml'

    # Open file containing bgc data
    with open(model_file_path, 'r') as f:
        model_info = yaml.full_load(f)
        base_element = model_info["base_element"]
        model = model_info["tracers"]
        reactions = model_info["reactions"]

    # base_element, reactions, tracers = import_bgc_model(model_file_path, physical)
    tracers = import_bgc_model(base_element, model, reactions, physical)

    concentration, sinking, tracer_map, tracer_type = create_function_inputs(physical["simulation"]["iters"],tracers)

    print('Model initialization complete.\n')

    return base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical