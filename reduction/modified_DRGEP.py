import copy
import logging
import numpy as np
from functions.bgc_rate_eqns import reduced_bgc_rate_eqns
import reduction.error_functions as error_functions
from reduction.pyMARS_DRGEP_functions import get_importance_coeffs
from scipy.integrate import solve_ivp
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List


def configure_reduction_dict(safe, target, tracer_map, tracers):

    tracer_names = []  # list of strings, names of all bgc tracers and their constituents

    safe_names = []     # list of strings, names of safe tracers and their constituents
    safe_indices = []   # list of ints, indices of safe tracers in concentration matrix
    
    target_names = []   # list of strings, names of target tracers and their constituents
    target_indices = [] # list of ints, indices of targets in concentration matrix

    for tracer,constituents in tracer_map.items():
        # Parse tracer dictionary
        for const in tracers[tracer].composition:
            # Add tracer constituent to list of species names
            tracer_names.append(f"{tracer}_{const}")

    if type(safe) is not type(None):
        for tracer,constituents in safe.items():
            # Add tracer to list of safe names ("all or nothing" scheme, all constituents for safe tracers are retained)
            # safe_names.append(tracer)
            # for i in range(len(tracer_map[tracer])):
            #     # Add index to list of safe indices
            #     safe_indices.append(tracer_map[tracer][i])

            for const in tracers[tracer].composition:
                if const in constituents:   
                    index = tracers[tracer].composition.index(const)
                    safe_names.append(f"{tracer}_{const}")
                    safe_indices.append(tracer_map[tracer][index])
    
    for tracer,constituents in target.items():
        for const in constituents:
            # Identify index of target constituent in tracer
            constituent_index = tracers[tracer].composition.index(const)
        
            # Extract associated index from tracer_map
            target_index = tracer_map[tracer][constituent_index]
        
            # Add index to list of target indices
            target_indices.append(target_index)

            # Add tracer constituent to list of target names
            target_names.append(f"{tracer}_{const}")

    # Create output dictionary
    dictionary = {}
    dictionary["tracer_names"] = tracer_names
    dictionary["safe_names"] = safe_names
    dictionary["safe_indices"] = safe_indices
    dictionary["target_names"] = target_names
    dictionary["target_indices"] = target_indices

    return dictionary


def slice_concentration(conc, physical):
    
    # Split depth to create slices
    depth_indices = np.array_split(np.arange(conc.shape[1]),6)
    indices = [split[0] for split in depth_indices]
    indices.append(depth_indices[-1][-1])

    # Create matrix for slices of concentration
    sliced_conc = conc[:,indices]

    # Create new dictionary for physical variables
    # sliced_physical = {
    #     "configuration": physical["simulation"]["configuration"],
    #     "dt": physical["simulation"]["dt"],
    #     "forcing": physical["environment"]["forcing"],
    #     "forcing_data": physical["environment"]["forcing_data"],
    #     "latitude": physical["environment"]["latitude"],
    #     "light_attenuation_water": physical["environment"]["light_attenuation_water"],
    #     "column_depth": physical["water_column"]["column_depth"],
    #     "num_boxes": len(indices),
    #     "z": physical["vertical_grid"]["z"][indices],
    #     "dz": physical["vertical_grid"]["dz"][indices],
    # }
    sliced_physical = {
        "configuration": physical["configuration"],
        "dt": physical["dt"],
        "forcing": physical["forcing"],
        "forcing_data": physical["forcing_data"],
        "latitude": physical["latitude"],
        "light_attenuation_water": physical["light_attenuation_water"],
        "column_depth": physical["column_depth"],
        "num_boxes": len(indices),
        "z": physical["z"][indices],
        "dz": physical["dz"][indices],
    }

    return sliced_conc, sliced_physical


def group_overall_interaction_coeffs(overall_interaction_coeffs, tracer_map, tracer_names):
    """ groups overall interaction coefficients by tracer
        coefficient for each tracer is the maximum value among its chemical constituents
    """

    for tracer,indices in tracer_map.items():
        # Determine maximum overall interaction coefficient among all chemical constituents in tracer
        max_coeff = 0.
        for index in indices:
            if overall_interaction_coeffs[tracer_names[index]] > max_coeff:     max_coeff = overall_interaction_coeffs[tracer_names[index]]

        # Apply maximum
        for index in indices:   overall_interaction_coeffs[tracer_names[index]] = max_coeff
        # overall_interaction_coeffs[indices] = overall_interaction_coeffs[indices].max()

    return overall_interaction_coeffs


def reduced_model_configuration(tracer_names, tracers_removed, old_conc, old_tracers):
    """ creates concentration matrix, tracer map, and tracer dictionary for reduced model
        removes tracers and associated reactions of tracers listed in "species_removed"
        reduces concentration matrix and number of functions called in simulation
    """

    # Reassign concentrations of removed tracers and split tracer names/constituents
    removed_tracer_names = []
    for name in tracers_removed:
        # Extract tracer name from list of removed species
        tracer,constituent = name.split("_")
        removed_tracer_names.append(tracer)

    # Remove repeated names
    removed_tracer_names = list(dict.fromkeys(removed_tracer_names))
    
    # Determine indices for removal
    indices_to_retain = []
    for index in range(0,len(tracer_names)):
        if tracer_names[index] not in tracers_removed: indices_to_retain.append(index)

    # Create new tracer dictionary
    new_tracers = copy.copy(old_tracers)
    for tracer in removed_tracer_names:
        # Remove from tracer dictionary
        if tracer in new_tracers:   del new_tracers[tracer]

    for tracer in list(new_tracers):
        for i in range(len(new_tracers[tracer].reactions)-1, -1, -1):
            # Create lists of consumed and produced keys
            if new_tracers[tracer].reactions[i]["consumed"] is not None:    consumed = list(new_tracers[tracer].reactions[i]["consumed"])
            else:   consumed = []   # Empty list if no consumed tracers
            if new_tracers[tracer].reactions[i]["produced"] is not None:    produced = list(new_tracers[tracer].reactions[i]["produced"])
            else:   produced = []   # Empty list if no produced tracers

            # Delete removed tracers from consumed list
            if any(name in consumed for name in removed_tracer_names):
                for name in consumed:
                    if name in removed_tracer_names:    new_tracers[tracer].reactions[i]["consumed"].pop(name, None)

            # Delete removed tracers from produced list
            if any(name in produced for name in removed_tracer_names):
                for name in produced:
                    if name in removed_tracer_names:    new_tracers[tracer].reactions[i]["produced"].pop(name, None)

            # Remove reactions from tracer reaction lists
            if type(new_tracers[tracer].reactions[i]["consumed"]) is not type(None) and type(new_tracers[tracer].reactions[i]["produced"]) is not type(None):
                if len(new_tracers[tracer].reactions[i]["consumed"]) == 0 and len(new_tracers[tracer].reactions[i]["produced"]) == 0:   
                    new_tracers[tracer].reactions.pop(i)
            
                elif len(new_tracers[tracer].reactions[i]["consumed"]) == 0 and len(new_tracers[tracer].reactions[i]["produced"]) != 0:
                    if new_tracers[tracer].reactions[i]["type"] not in ["chlorophyll_synthesis","co2_flux","gross_primary_production","reaeration"]:   
                        new_tracers[tracer].reactions.pop(i)
            
                elif len(new_tracers[tracer].reactions[i]["consumed"]) != 0 and len(new_tracers[tracer].reactions[i]["produced"]) == 0:
                    if new_tracers[tracer].reactions[i]["type"] not in ["denitrification","nitrification","photosynthesis","respiration"]: 
                        new_tracers[tracer].reactions.pop(i)

    # Set empty dicts to 'None'
    for tracer in list(new_tracers):
        for i in range(len(new_tracers[tracer].reactions)-1, -1, -1):
            if isinstance(new_tracers[tracer].reactions[i]["consumed"],dict) and len(new_tracers[tracer].reactions[i]["consumed"]) == 0:
                new_tracers[tracer].reactions[i]["consumed"] = None
            if isinstance(new_tracers[tracer].reactions[i]["produced"],dict) and len(new_tracers[tracer].reactions[i]["produced"]) == 0:
                new_tracers[tracer].reactions[i]["produced"] = None

    # Create new tracer map
    new_tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
    new_tracer_type = []   # used in vertical diffusivity calculations
    index = 0
    for trac in new_tracers:
        num_constituents = len(new_tracers[trac].composition)   # number of constituents in tracer
        
        lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
        for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
        new_tracer_map[trac] = lst  # identify tracer constituents with their own index
        
        for i in range(num_constituents):
            # add tracer type to list
            if new_tracers[trac].type == "detritus":    new_tracer_type.append(new_tracers[trac].form)     # need to distinguish particulate/dissolved form
            else:   new_tracer_type.append(new_tracers[trac].type)     # just the type
        
            # add tracer type to list
            index += 1  # update index

    # Create new concentration matrix
    new_conc = np.zeros((len(indices_to_retain),old_conc.shape[1]))
    for i in range(0,len(indices_to_retain)):
        new_conc[i,:] = old_conc[indices_to_retain[i],:]
    
    return new_tracers, new_conc


# def reduced_model_configuration(tracer_names, tracers_removed, old_conc, old_tracer_map, old_tracers, reassign):
#     """ creates concentration matrix, tracer map, and tracer dictionary for reduced model
#         removes tracers and associated reactions of tracers listed in "species_removed"
#         reduces concentration matrix and number of functions called in simulation
#     """

#     # Reassign concentrations of removed tracers and split tracer names/constituents
#     removed_tracer_names = []
#     for name in tracers_removed:
#         # Extract tracer name from list of removed species
#         tracer,constituent = name.split("_")
#         removed_tracer_names.append(tracer)

#         # Determine tracer type, use form if type is detritus
#         if old_tracers[tracer].type in ["bacteria","phytoplankton","zooplankton","detritus"]:
#             if old_tracers[tracer].type == "detrirus":  tracer_type = old_tracers[tracer].form
#             else:   tracer_type = old_tracers[tracer].type

#             # Parse reassign dictionary
#             if type(reassign) is not type(None):
#                 for key,val in reassign.items():
#                     if tracer_type == key and tracer in val:
#                         # Extract compositions and associated indices in concentration matrix
#                         original_composition = old_tracers[tracer].composition      # Composition of tracer being removed
#                         original_indices = old_tracer_map[tracer]                   # Indices of removed tracer in concentration matrix

#                         # Index of first tracer that isn't the removed tracer
#                         # ex: val = ["phyto1","phyto2","phyto3","phyto4"] --> if tracer = "phyto1" then index = 1, if tracer != "phyto1" then index = 0
#                         index = next(i for i,name in enumerate(val) if name != tracer)  
#                         reassign_composition = old_tracers[val[index]].composition  # Composition of tracer getting the reassigned concentration
#                         reassign_indices = old_tracer_map[val[index]]               # Indices of reassigned tracer in concentration matrix

#                         # Reassign concentrations by chemical constituent
#                         for const in original_composition:
#                             if const in reassign_composition:
#                                 # Get constituent index of both tracers
#                                 oc_const_index = original_indices[original_composition.index(const)]
#                                 rc_const_index = reassign_indices[reassign_composition.index(const)]

#                                 # Reassign concentration
#                                 old_conc[rc_const_index] += old_conc[oc_const_index]

#                         # Delete removed tracer from reassign dictionary
#                         reassign[key].remove(tracer)
#                         # val.remove(tracer)

#     # Remove repeated names
#     removed_tracer_names = list(dict.fromkeys(removed_tracer_names))
    
#     # Determine indices for removal
#     indices_to_retain = []
#     for index in range(0,len(tracer_names)):
#         if tracer_names[index] not in tracers_removed: indices_to_retain.append(index)

#     # Create new tracer dictionary
#     new_tracers = copy.copy(old_tracers)
#     for tracer in removed_tracer_names:
#         # Remove from tracer dictionary
#         if tracer in new_tracers:   del new_tracers[tracer]

#     for tracer in list(new_tracers):
#         for i in range(len(new_tracers[tracer].reactions)-1, -1, -1):
#             # Create lists of consumed and produced keys
#             if new_tracers[tracer].reactions[i]["consumed"] is not None:    consumed = list(new_tracers[tracer].reactions[i]["consumed"])
#             else:   consumed = []   # Empty list if no consumed tracers
#             if new_tracers[tracer].reactions[i]["produced"] is not None:    produced = list(new_tracers[tracer].reactions[i]["produced"])
#             else:   produced = []   # Empty list if no produced tracers

#             # Delete removed tracers from consumed list
#             if any(name in consumed for name in removed_tracer_names):
#                 for name in consumed:
#                     if name in removed_tracer_names:    new_tracers[tracer].reactions[i]["consumed"].pop(name, None)

#             # Delete removed tracers from produced list
#             if any(name in produced for name in removed_tracer_names):
#                 for name in produced:
#                     if name in removed_tracer_names:    new_tracers[tracer].reactions[i]["produced"].pop(name, None)

#             # Remove reactions from tracer reaction lists
#             if type(new_tracers[tracer].reactions[i]["consumed"]) is not type(None) and type(new_tracers[tracer].reactions[i]["produced"]) is not type(None):
#                 if len(new_tracers[tracer].reactions[i]["consumed"]) == 0 and len(new_tracers[tracer].reactions[i]["produced"]) == 0:   
#                     new_tracers[tracer].reactions.pop(i)
            
#                 elif len(new_tracers[tracer].reactions[i]["consumed"]) == 0 and len(new_tracers[tracer].reactions[i]["produced"]) != 0:
#                     if new_tracers[tracer].reactions[i]["type"] not in ["chlorophyll_synthesis","co2_flux","gross_primary_production","reaeration"]:   
#                         new_tracers[tracer].reactions.pop(i)
            
#                 elif len(new_tracers[tracer].reactions[i]["consumed"]) != 0 and len(new_tracers[tracer].reactions[i]["produced"]) == 0:
#                     if new_tracers[tracer].reactions[i]["type"] not in ["denitrification","nitrification","photosynthesis","respiration"]: 
#                         new_tracers[tracer].reactions.pop(i)

#     # Set empty dicts to 'None'
#     for tracer in list(new_tracers):
#         for i in range(len(new_tracers[tracer].reactions)-1, -1, -1):
#             if isinstance(new_tracers[tracer].reactions[i]["consumed"],dict) and len(new_tracers[tracer].reactions[i]["consumed"]) == 0:
#                 new_tracers[tracer].reactions[i]["consumed"] = None
#             if isinstance(new_tracers[tracer].reactions[i]["produced"],dict) and len(new_tracers[tracer].reactions[i]["produced"]) == 0:
#                 new_tracers[tracer].reactions[i]["produced"] = None

#     # Create new tracer map
#     new_tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
#     new_tracer_type = []   # used in vertical diffusivity calculations
#     index = 0
#     for trac in new_tracers:
#         num_constituents = len(new_tracers[trac].composition)   # number of constituents in tracer
        
#         lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
#         for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
#         new_tracer_map[trac] = lst  # identify tracer constituents with their own index
        
#         for i in range(num_constituents):
#             # add tracer type to list
#             if new_tracers[trac].type == "detritus":    new_tracer_type.append(new_tracers[trac].form)     # need to distinguish particulate/dissolved form
#             else:   new_tracer_type.append(new_tracers[trac].type)     # just the type
        
#             # add tracer type to list
#             index += 1  # update index

#     # Create new concentration matrix
#     new_conc = np.zeros((len(indices_to_retain),old_conc.shape[1]))
#     for i in range(0,len(indices_to_retain)):
#         new_conc[i,:] = old_conc[indices_to_retain[i],:]
    

#     # return new_conc, new_tracers, new_tracer_map, new_tracer_type, reassign
#     return new_tracers, new_conc


def calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t, base_element, physical, tracers, tracer_map, tracer_type):
    """ Calculates the percent difference between the new rate eqn and the original.
        The new rate is based on turning one species 'off' """

    c_original = copy.copy(conc)
    c_new = copy.copy(conc)
    num_tracers = conc.shape[0]
    num_boxes = conc.shape[1]

    indices_to_retain = np.arange(num_tracers)
    removed_tracer_names = []

    # calculate original rate values
    # dc_dt_og = rate_eqn_fcn(t, base_element, conc, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names, True)
    dc_dt_og = rate_eqn_fcn(t, base_element, conc, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names)

    percent_error_matrix = np.zeros([num_tracers, num_tracers, num_boxes])

    for j in range(num_tracers):
        c_new[j,:] = 0.     # zero out concentration of current tracer
        
        # dc_dt_new = rate_eqn_fcn(t, base_element, c_new, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names, True)
        dc_dt_new = rate_eqn_fcn(t, base_element, c_new, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names)
        # new = dc_dt_new
        # og = dc_dt_og
        new = dc_dt_new.reshape((num_tracers,num_boxes))
        og = dc_dt_og.reshape((num_tracers,num_boxes))
        percent_error = calc_percent_error(new, og, num_tracers, num_boxes)
        percent_error_matrix[:,j,:] = percent_error
        c_new[j,:] = c_original[j,:]    # restore original concentration

            
    percent_error = np.amax(percent_error_matrix,axis=2)

    # Find the maximum value along each row (used for normalization)
    row_max = np.amax(percent_error, axis=1)

    # Normalize percent_error_matrix by row max which is new_dic_matrix
    dic_matrix = np.zeros([num_tracers, num_tracers])
    for i in range(num_tracers):
        if row_max[i] == 0:
            dic_matrix[i,:] = 0.0
        else:
            dic_matrix[i,:] = percent_error[i,:]/row_max[i]

    # Set diagonals to zero to avoid self_directing graph edges
    np.fill_diagonal(dic_matrix, 0.0) 

    return dic_matrix


def calc_percent_error(new_matrix,old_matrix,num_tracers,num_boxes):
    """ Calculates the percent error between two matricies.
    This is used for the calculating the New Method's direct interaction coeffs
    """
    percent_error = np.zeros_like(new_matrix)
    if new_matrix.shape == old_matrix.shape:
        for i in range(0,num_tracers):
            for k in range(0,num_boxes):
                if new_matrix[i,k] != old_matrix[i,k]:
                    percent_error[i,k] = 100*(abs(new_matrix[i,k] - old_matrix[i,k])/abs(old_matrix[i,k]))
                if np.isnan(percent_error[i,k]):
                    percent_error[i,k] = 0


    return percent_error


def modified_DRGEP(conc, reduction, base_element, physical, tracer_map, tracer_type, tracers):
    """
    Description: Apply Modififed DRGEP reduction strategy to BGC model.
                 Target and safe species for each error function listed below.
    """
    
    # Information for reduction
    scenario = reduction["title"]
    error_limit = reduction["error_limit"]

    # Build dictionary of original model reduction information
    og_reduction_info = configure_reduction_dict(reduction["safe"], reduction["target"], tracer_map, tracers)
    # og_reduction_info["reassign"] = reduction["reassign"]
    og_reduction_info["time_period"] = reduction["time_period"]
    tracer_names = og_reduction_info["tracer_names"]
    target_names = og_reduction_info["target_names"]
    safe_names = og_reduction_info["safe_names"]
    if reduction["mode"] == "average":
        error_function = error_functions.average
    elif reduction["mode"] == "peak":
        error_function = error_functions.peak
    elif reduction["mode"] == "time_of_peak":
        error_function = error_functions.time_of_peak
    rate_eqn_fcn = reduced_bgc_rate_eqns

    # Time span for integration
    # t_span = [0,86400*365*10]
    t_span = [0,86400*365*3]

    # Time at which the DIC values are obtained
    t = 86400*365*7

    # Log input data to file
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename='output.log', level=logging.INFO)
    logging.info(100 * '-')
    logging.info('Scenario: {}'.format(scenario))
    logging.info('Error function: {}'.format(error_function))
    logging.info('Error limit: {}'.format(error_limit))
    logging.info('Target species: {}'.format(target_names))
    logging.info('Retained species: {}'.format(safe_names))

    # dic_physical = {
    #     "configuration": physical["simulation"]["configuration"],
    #     "dt": physical["simulation"]["dt"],
    #     "forcing": physical["environment"]["forcing"],
    #     "forcing_data": physical["environment"]["forcing_data"],
    #     "latitude": physical["environment"]["latitude"],
    #     "light_attenuation_water": physical["environment"]["light_attenuation_water"],
    #     "column_depth": physical["water_column"]["column_depth"],
    #     "num_boxes": physical["water_column"]["num_boxes"],
    #     "z": physical["vertical_grid"]["z"],
    #     "dz": physical["vertical_grid"]["dz"],
    # }

    # Get direct interaction coefficients
    # dic_matrix = calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t, base_element, dic_physical, tracers, tracer_map, tracer_type)
    dic_matrix = calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t, base_element, physical, tracers, tracer_map, tracer_type)

    # Get overall interaction coefficients
    overall_interaction_coeffs = get_importance_coeffs(tracer_names, target_names, [dic_matrix])

    # Group overall interaction coefficients
    overall_interaction_coeffs = group_overall_interaction_coeffs(overall_interaction_coeffs, tracer_map, tracer_names)

    # Slice concentration matrix
    if conc.shape[1] > 7:  # separate into 6 slices (top, bottom, and 4 interior slices - runs faster than inputting full concentration matrix)
        sliced_conc, sliced_physical = slice_concentration(conc, physical)
    else:
        sliced_conc = copy.copy(conc)
        sliced_physical = copy.copy(physical)
        # sliced_physical = {
        #     "configuration": physical["simulation"]["configuration"],
        #     "dt": physical["simulation"]["dt"],
        #     "forcing": physical["environment"]["forcing"],
        #     "forcing_data": physical["environment"]["forcing_data"],
        #     "latitude": physical["environment"]["latitude"],
        #     "light_attenuation_water": physical["environment"]["light_attenuation_water"],
        #     "column_depth": physical["water_column"]["column_depth"],
        #     "num_boxes": physical["water_column"]["num_boxes"],
        #     "z": physical["vertical_grid"]["z"],
        #     "dz": physical["vertical_grid"]["dz"],
        # }
    
    # Run new method
    reduction_data = run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_function, og_reduction_info, t_span, sliced_conc, base_element, sliced_physical, tracers, tracer_map, tracer_type)

    return reduction_data, error_limit


def reduce_modified_DRGEP(error_function, last_error, tracer_names, target_names, safe_names, threshold, overall_interaction_coeffs, last_num_tracers, solution_full_model, solution_reduced_model, t_span, time_period, base_element, c0, tracers, tracer_map, tracer_type, physical):
    """ calculates the number of species and error for a given threshold value
    """

    # Find species to remove using cutoff threshold
    tracers_removed = []
    for tracer,coeff in overall_interaction_coeffs.items():
        if coeff < threshold and tracer not in safe_names:   tracers_removed.append(tracer)

    # Reassign concentrations of removed tracers and split tracer names/constituents
    removed_tracer_names = []
    for name in tracers_removed:
        # Extract tracer name from list of removed species
        tracer,constituent = name.split("_")
        removed_tracer_names.append(tracer)
    
    # Remove repeated names
    removed_tracer_names = list(dict.fromkeys(removed_tracer_names))
    
    # Determine indices for removal
    indices_to_retain = []
    for index in range(0,len(tracer_names)):
        if tracer_names[index] not in tracers_removed: indices_to_retain.append(index)
    
    for index in range(0,c0.shape[0]):
        if index not in indices_to_retain:  c0[index,...] = 0.

    # Count how many species are remaining
    num_tracers = len(tracer_names) - len(tracers_removed)

    if num_tracers == last_num_tracers:
        error = last_error
        solution_reduced_model = solution_reduced_model
    else:
       # For the reduced model, call function to get error
        error, solution_reduced_model = error_function(t_span, time_period, base_element, physical, tracer_names, target_names, solution_full_model, c0, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names)

    return error, num_tracers, tracers_removed, solution_reduced_model, c0, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names


def run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_function, og_reduction_info, t_span, c0_full, base_element, physical, tracers_full, tracer_map_full, tracer_type_full):
    
    """ Iterates through different threshold values to find a reduced model that meets error criteria
    """

    tracer_names = og_reduction_info["tracer_names"]
    target_names = og_reduction_info["target_names"]
    safe_names = og_reduction_info["safe_names"]

    assert target_names, 'Need to specify at least one target species.'

    # begin reduction iterations
    logging.info('Beginning reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    # make lists to store data
    threshold_data = []
    num_tracers_data = []
    error_data = []
    tracers_removed_data = []
    solution_reduced_models = []
    output_data = {}

    # start with detailed (starting) model
    num_tracers = c0_full.shape[0]
    indices_to_retain = np.arange(num_tracers)
    removed_tracer_names = []
    # solution_full_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers, physical, tracers_full, tracer_map_full, tracer_type_full, indices_to_retain, removed_tracer_names, False), t_span, c0_full.ravel(), method='RK23')#, max_step=physical["dt"])
    solution_full_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers, physical, tracers_full, tracer_map_full, tracer_type_full, indices_to_retain, removed_tracer_names), t_span, c0_full.ravel(), method='RK23')#, max_step=physical["dt"])

    first = True
    error_current = 0.0
    threshold = 5e-2#1e-3#3e-2
    threshold_increment = 3e-3#1e-3 #1e-5
    threshold_multiplier = 2
    while error_current <= error_limit:
        if first: 
            solution_reduced_model = solution_full_model
            c0 = c0_full
            tracers = tracers_full
            tracer_map = tracer_map_full
            tracer_type = tracer_type_full
            indices_to_retain = []
            removed_tracer_names = []
        error_current, num_tracers, species_removed, solution_reduced_model, c0, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names = reduce_modified_DRGEP(error_function, error_current, tracer_names, target_names, safe_names, threshold, overall_interaction_coeffs, num_tracers, 
                                                                                                                                                                         solution_full_model, solution_reduced_model, t_span, og_reduction_info["time_period"], base_element, c0, tracers, tracer_map, tracer_type, physical)

        # reduce threshold if past error limit on first iteration
        if first and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-10:
                raise SystemExit(
                    'Threshold value dropped below 1e-10 without producing viable reduced model'
                    )
            logging.info('Threshold value too high, reducing by factor of 10')
            continue

        logging.info(f'{threshold:^9.2e} | {num_tracers:^17} | {error_current:^.5f}')

        # store data
        threshold_data.append(threshold)
        num_tracers_data.append(num_tracers)
        error_data.append(error_current)
        tracers_removed_data.append(species_removed)
        solution_reduced_models.append(solution_reduced_model)

        threshold += threshold_increment
        # threshold *= threshold_multiplier
        first = False

        # Stop reduction process if num species reaches one
        if num_tracers <= 1:
            break

        # Stop interating if threshold exceeds value
        # if threshold >= 2e-3:
        # if threshold >= 3e-2:
        if threshold >= 1:  # maximum threshold for oxy_3 error function   
            break

    if error_current > error_limit:
        threshold -= (2 * threshold_increment)
        error_current, num_tracers, species_removed, solution_reduced_model, c0, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names = reduce_modified_DRGEP(error_function, error_current, tracer_names, target_names, safe_names, threshold, overall_interaction_coeffs, num_tracers, 
                                                                                                                                                                         solution_full_model, solution_reduced_model, t_span, og_reduction_info["time_period"], base_element, c0, tracers, tracer_map, tracer_type, physical)
    
    if error_data[-1] > error_limit:
        model = -2
    else:
        model = -1 
        
    # Make dictionary of species remaining in reduced model
    remaining_species = {}
    for index, species in enumerate(tracer_names):
        if species not in tracers_removed_data[model]:
            remaining_species[species] = index
          
    logging.info(45 * '-')
    logging.info('New method reduction complete.')
    logging.info('Remaining species: {}'.format(remaining_species))
    logging.info('Removed species: {}'.format(tracers_removed_data[model]))
    logging.info(100 * '-')

    # Store all output data to a dictionary
    output_data['threshold_data'] = threshold_data
    output_data['num_tracers_data'] = num_tracers_data
    output_data['error_data'] = error_data
    output_data['tracers_removed_data'] = tracers_removed_data
    output_data['solution_full_model'] = solution_full_model
    output_data['solution_reduced_models'] = solution_reduced_models

    output_data['tracer_names'] = tracer_names
    output_data['tracers_reduced'] = tracers
    output_data['tracer_map_reduced'] = tracer_map
    output_data['tracer_type_reduced'] = tracer_type
    output_data['conc_reduced'] = c0

    return output_data
