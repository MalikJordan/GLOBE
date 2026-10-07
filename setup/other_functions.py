import copy
import numpy as np
import yaml
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List

class BlankNoneDumper(yaml.SafeDumper):
    pass


def represent_none(dumper, value):
    return dumper.represent_scalar(
        "tag:yaml.org,2002:null",
        ""
    )


def create_function_inputs(iters, tracers):
    """
    Definition: Takes tracer dictionary and creates lists, arrays, or typed.Dicts for numba calculations

    :return: concentration (array), sinking velocities (array), tracer map (typed.Dict), tracer types (list)
    """

    # Create list of concentrations
    initial_concentration = []

    # Create typed.Dict of tracer indices in concentration
    tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
    
    # Create list of trcaer types
    tracer_type = []   # used in vertical diffusivity calculations

    # Create list of sinking velocities for each tracer
    sinking = []

    index = 0   # counting number for tracer indices
    for trac in tracers:
        num_constituents = len(tracers[trac].composition)   # number of constituents in tracer

        lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
        for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
        tracer_map[trac] = lst  # identify tracer constituents with their own index

        for i in range(num_constituents):
            # add concentration to matrix
            initial_concentration.append(tracers[trac].initial_conc[i,...])    # add concentration to matrix

            # add tracer type to list
            if tracers[trac].type == "detritus":    tracer_type.append(tracers[trac].form)     # need to distinguish particulate/dissolved form
            else:   tracer_type.append(tracers[trac].type)     # just the type

            # add sinking velocity to list
            if hasattr(tracers[trac],"sinking_velocity"):   sinking.append(tracers[trac].sinking_velocity)
            else:   sinking.append(np.zeros(tracers[trac].initial_conc.shape[1]))

            # add tracer type to list
            index += 1  # update index

    initial_concentration = np.array(initial_concentration,dtype=np.float64)    # convert concentration from list to array
    concentration = np.zeros((initial_concentration.shape[0],initial_concentration.shape[1],iters),dtype=np.float64)
    concentration[:,:,0] = initial_concentration.copy()
    sinking = np.array(sinking,dtype=np.float64)    # convert sinking from list to array

    return concentration, sinking, tracer_map, tracer_type


def clean(obj, removed_tracer_names, parameter_fields, parent_key=None):
    """
    Definition:: Recursive function to clean parameter space of removed tracers
    """
    if isinstance(obj, dict):
        cleaned = {}    # Initialize cleaned parameter dictionary

        for key, value in obj.items():
            if parent_key in parameter_fields and key in removed_tracer_names:  continue
            cleaned_value = clean(value, removed_tracer_names, parameter_fields, key)
            if cleaned_value is None:   continue
            if isinstance(cleaned_value, (dict, list)) and not cleaned_value:   continue

            cleaned[key] = cleaned_value

        return cleaned

    elif isinstance(obj, list):
        cleaned = []    # Initialize cleaned parameter list

        for value in obj:
            if parent_key in parameter_fields and isinstance(value, str) and value in removed_tracer_names: continue
            cleaned_value = clean(value, removed_tracer_names, parameter_fields, parent_key)
            if cleaned_value is not None:   cleaned.append(cleaned_value)

        return cleaned

    else:
        if parent_key in parameter_fields and isinstance(obj, str) and obj in removed_tracer_names:
            return None

        return obj
    
        
def trim_yaml(file_path, tracers_removed, reduced_tracers):
    """
    Definition:: Takes original input yaml file and trims removed tracers and associated parameters
    
    :return: Updated yaml file
    """
    # ----------------------------------------------------------------------------------------------------
    # Open file containing bgc data
    # ----------------------------------------------------------------------------------------------------
    with open(file_path, 'r') as f:
        model_info = yaml.full_load(f)
        base_element = model_info["base_element"]
        old_tracers = model_info["tracers"]

    # ----------------------------------------------------------------------------------------------------
    # Remove tracers
    # ----------------------------------------------------------------------------------------------------
    # Split names/constituents for removed tracers
    removed_tracer_names = []
    for name in tracers_removed:
        # Extract tracer name from list of removed species
        tracer,constituent = name.split("_")
        removed_tracer_names.append(tracer)

    # Remove repeated names
    removed_tracer_names = list(dict.fromkeys(removed_tracer_names))
    
    # Delete removed tracers from dictionary
    new_tracers = copy.copy(old_tracers)    # Create new tracer dictionary
    for tracer in removed_tracer_names:
        if tracer in new_tracers:   del new_tracers[tracer]     # Remove from tracer dictionary

    # ----------------------------------------------------------------------------------------------------
    # Rewrite reactions list
    # ----------------------------------------------------------------------------------------------------
    reactions = []  # Initialize empty list
    for tracer in reduced_tracers:
        for reac in reduced_tracers[tracer].reactions:    reactions.append(reac)
    
    # ----------------------------------------------------------------------------------------------------
    # Clean old YAML input file
    # ----------------------------------------------------------------------------------------------------
    # Define parameter fields for necessary parsing
    parameter_fields = {"grazing_preferences", "om_partition", "uptake", "nutrient_limitation", "substrates", "potential_rich", "potential_poor"}
    
    # Clean parameters
    for tracer in new_tracers.values():
        if "parameters" in tracer:  tracer["parameters"] = clean(tracer["parameters"], removed_tracer_names, parameter_fields)

    # ----------------------------------------------------------------------------------------------------
    # Create new YAML input file
    # ----------------------------------------------------------------------------------------------------
    # Create new configuration dictionary
    configuration = {
        "base_element": base_element,
        "tracers": new_tracers,
        "reactions": reactions
    }

    # Keep "None" as empty instead of writing "null"
    BlankNoneDumper.add_representer(type(None), represent_none)

    # Write configuration to new YAML input file
    with open("reduced_model.yaml", "w") as file:
        yaml.dump(configuration, file, Dumper=BlankNoneDumper, sort_keys=False, default_flow_style=False)
