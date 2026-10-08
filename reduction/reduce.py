import copy
import os
import yaml
from setup.other_functions import create_function_inputs, trim_yaml
from reduction.modified_DRGEP import modified_DRGEP, reduced_model_configuration


def reduce_bgc_model(model_file_path, physical, base_element, concentration, tracer_map, tracer_type, tracers):
    """
    Definition:: handles reduction scheme and rewrites YAML input file with reduced model configuration
    """
    with open(os.getcwd() + '/reduction.yaml', 'r') as f:
        reduction_scheme = yaml.full_load(f)

    # Create dictionary of physical variables
    reduction_physical = {
        "configuration": physical["simulation"]["configuration"],
        "dt": physical["simulation"]["dt"],
        "num_boxes": physical["water_column"]["num_boxes"],
        "column_depth": physical["water_column"]["column_depth"],
        "latitude": physical["environment"]["latitude"],
        "light_attenuation_water": physical["environment"]["light_attenuation_water"],
    }

    # Add environmental forcing
    if physical["environment"]["forcing"] == "pom1d":
        reduction_physical["forcing"] = "seasonal"
        reduction_physical["forcing_data"] = {
            "summer_mld": 10.0,
            "winter_mld": 40.0,
            "summer_salt": 36.5,
            "winter_salt": 37.,
            "summer_sun": 300.0,
            "winter_sun": 20.0,
            "summer_temp": 30.,
            "winter_temp": 10.,
            "temp_excursion": 1.,
            "summer_wind": 2.,
            "winter_wind": 6.
        }
    else:
        reduction_physical["forcing"] = copy.copy(physical["environment"]["forcing"])
        reduction_physical["forcing_data"] = copy.copy(physical["environment"]["forcing_data"])

    if physical["simulation"]["configuration"] == "1d":
        reduction_physical["z"] = copy.copy(physical["vertical_grid"]["z"])[:-1]
        reduction_physical["dz"] = copy.copy(physical["vertical_grid"]["dz"])[:-1]
    else:
        reduction_physical["z"] = copy.copy(physical["vertical_grid"]["z"])
        reduction_physical["dz"] = copy.copy(physical["vertical_grid"]["dz"])

    # Reduce model
    reduction, error_limit = modified_DRGEP(concentration[...,0].copy(), reduction_scheme, base_element, reduction_physical, tracer_map, tracer_type, tracers)

    error_data = reduction['error_data'][-2:]
    if error_data[-1] > error_limit:    model = -2
    else:   model = -1

    # Extract data
    tracers_removed = reduction["tracers_removed_data"][model]
    tracer_names = reduction["tracer_names"]
    
    # BFM40
    # tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'ta_eq']
    # BFM38
    # tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'dom2_c', 'dom3_c', 'ta_eq']
    # BFM37
    # tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'dom2_c', 'dom3_c', 'co2_c', 'ta_eq']
    # BFM24
    # tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'phyto1_c', 'phyto1_n', 'phyto1_p', 'phyto1_chl', 'phyto1_si', 'phyto3_c', 'phyto3_n', 'phyto3_p', 'phyto3_chl', 'phyto4_c', 'phyto4_n', 'phyto4_p', 'phyto4_chl', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'microzoo2_c', 'microzoo2_n', 'microzoo2_p', 'ta_eq']
    # BFM23
    tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'phyto1_c', 'phyto1_n', 'phyto1_p', 'phyto1_chl', 'phyto1_si', 'phyto3_c', 'phyto3_n', 'phyto3_p', 'phyto3_chl', 'phyto4_c', 'phyto4_n', 'phyto4_p', 'phyto4_chl', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'microzoo2_c', 'microzoo2_n', 'microzoo2_p', 'dom3_c', 'ta_eq']
    # BFM22
    # tracers_removed = ['n2_n', 'sio4_si', 'hs_s', 'phyto1_c', 'phyto1_n', 'phyto1_p', 'phyto1_chl', 'phyto1_si', 'phyto3_c', 'phyto3_n', 'phyto3_p', 'phyto3_chl', 'phyto4_c', 'phyto4_n', 'phyto4_p', 'phyto4_chl', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'microzoo2_c', 'microzoo2_n', 'microzoo2_p', 'dom2_c', 'dom3_c', 'ta_eq']

    # Reconfigure model
    tracer_names = ['o2_o', 'po4_p', 'no3_n', 'nh4_n', 'n2_n', 'sio4_si', 'hs_s', 'bac1_c', 'bac1_n', 'bac1_p', 'phyto1_c', 'phyto1_n', 'phyto1_p', 'phyto1_chl', 'phyto1_si', 'phyto2_c', 'phyto2_n', 'phyto2_p', 'phyto2_chl', 'phyto3_c', 'phyto3_n', 'phyto3_p', 'phyto3_chl', 'phyto4_c', 'phyto4_n', 'phyto4_p', 'phyto4_chl', 'mesozoo1_c', 'mesozoo1_n', 'mesozoo1_p', 'mesozoo2_c', 'mesozoo2_n', 'mesozoo2_p', 'microzoo1_c', 'microzoo1_n', 'microzoo1_p', 'microzoo2_c', 'microzoo2_n', 'microzoo2_p', 'dom1_c', 'dom1_n', 'dom1_p', 'dom2_c', 'dom3_c', 'pom1_c', 'pom1_n', 'pom1_p', 'pom1_si', 'co2_c', 'ta_eq']
    reduced_conc, reduced_tracers, reduced_tracer_map, reduced_tracer_type = reduced_model_configuration(tracer_names, tracers_removed, concentration[...,0].copy(), tracers)
    
    trim_yaml(model_file_path, tracers_removed, reduced_tracers)