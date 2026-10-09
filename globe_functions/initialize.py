import os
from setup.other_functions import create_function_inputs
from setup.initialize import import_bgc_model, import_physical_model

def initialize_model():
    """
    Definition:: handles the initialization of bgc model and physical/environmental parameters
    """

    print('Initializing model...')
    # check_file = False
    # first_check = True
    # while not check_file:
    #     if first_check:
    #         file = input("\nEnter the name of the yaml input file. Press 'Q' to quit. \nEx: model_description.yaml \n\n")
    #     else:
    #         file = input("\nInput file not found. Re-enter the name of the yaml input file or press 'Q' to quit. \nEx: model_description.yaml \n\n")
    #     if file == 'Q' or file == 'q':
    #         sys.exit()
    #     check_file = os.path.exists(file)
    #     first_check = False

    # Import physical model
    # file = 'physical_bfm17_1d.yaml'
    # file = 'tests/bfm56/data/physical_bfm56.yaml'
    # physical_file_path = os.getcwd() + '/' + file
    physical_file_path = os.getcwd() + '/physical.yaml'
    physical = import_physical_model(physical_file_path)

    # file = 'bfm17_1d.yaml'
    # file = 'tests/bfm56/data/bfm56.yaml'
    # model_file_path = os.getcwd() + '/' + file
    model_file_path = os.getcwd() + '/model.yaml'
    base_element, reactions, tracers = import_bgc_model(model_file_path, physical)

    concentration, sinking, tracer_map, tracer_type = create_function_inputs(physical["simulation"]["iters"],tracers)

    print('Model initialization complete.\n')

    return base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical