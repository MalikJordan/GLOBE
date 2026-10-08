import numpy as np
import copy
from scipy.integrate import solve_ivp
from functions.bgc_rate_eqns import reduced_bgc_rate_eqns

def num_tracers(tracer_map):

    total = sum(len(v) for v in tracer_map.values())

    return total


def average(t_span, time_period, base_element, physical, tracer_names, target_names, solution_full_model, c0_reduced, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names):
    """ calculates error in the time of peak concentration of target tracer (sum of targets if multiple are provided)
    """
    num_tracers_full = len(tracer_names)
    num_tracers_reduced = num_tracers(tracer_map_reduced)
    
    num_boxes = c0_reduced.shape[1]     # number of boxes in water column
    
    # solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names, False), 
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names), 
                                 t_span, c0_reduced.ravel(), method='RK23')#, max_step=physical["dt"])
    
    # Determine start and end time for slicing
    t_min = 86400 * time_period["start"]
    t_max = 86400 * time_period["end"]
    
    # Find index associated with t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
    
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # Reshape concentration matrices
    concentration_full = solution_full_model.y.reshape(num_tracers_full,num_boxes,len(solution_full_model.t))               # shape = (num_tracers, num_boxes, time)
    concentration_reduced = solution_reduced_model.y.reshape(num_tracers_reduced,num_boxes,len(solution_reduced_model.t))   # shape = (num_tracers, num_boxes, time)
    
    # Find target indices
    target_indices_full = []
    target_indices_reduced = []
    
    for name in tracer_names:
        # Index in solution_full_model
        if name in target_names: target_indices_full.append(tracer_names.index(name))
    
        # Split tracer name and constituent
        tracer,constituent = name.split("_")
    
        # Index in solution_reduced_model
        if tracer in tracer_map_reduced:
            if name in target_names:
                # Index of target constituent in target tracer
                const_index = tracers_reduced[tracer].composition.index(constituent)
                target_indices_reduced.append(tracer_map_reduced[tracer][const_index])
    
    # Extract concentration over the time period
    target_quality_full = np.zeros((len(target_indices_full),num_boxes,(index_t_max_full-index_t_min_full)))                # shape = (num_targets, num_boxes, time)
    target_quality_reduced = np.zeros((len(target_indices_reduced),num_boxes,(index_t_max_reduced-index_t_min_reduced)))    # shape = (num_targets, num_boxes, time)
    
    for index_full in range(index_t_min_full,index_t_max_full):
        target_quality_full[...,(index_full - index_t_min_full)] = concentration_full[target_indices_full,:,index_full]
    for index_reduced in range(index_t_min_reduced,index_t_max_reduced):
        target_quality_reduced[...,(index_reduced - index_t_min_reduced)] = concentration_reduced[target_indices_reduced,:,index_reduced]

    # Collapse target quality matrix
    collapsed_full = np.sum(target_quality_full, axis=0, keepdims=True)         # shape = (1, num_boxes, time)
    collapsed_reduced = np.sum(target_quality_reduced, axis=0, keepdims=True)   # shape = (1, num_boxes, time)
    
    # Find the index of peak concentration
    average_full = np.mean(collapsed_full)
    average_reduced = np.mean(collapsed_reduced)
    
    # Compute the error with respect to full solution
    error = np.abs(100 * np.abs(average_full - average_reduced) / average_full)
        
    return error, solution_reduced_model


def peak(t_span, time_period, base_element, physical, tracer_names, target_names, solution_full_model, c0_reduced, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names):
    """ calculates error in the peak concentration of target tracer (sum of targets if multiple are provided)
    """
    num_tracers_full = len(tracer_names)
    num_tracers_reduced = num_tracers(tracer_map_reduced)

    num_boxes = c0_reduced.shape[1]     # number of boxes in water column

    # solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names, False), 
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names), 
                                 t_span, c0_reduced.ravel(), method='RK23')#, max_step=physical["dt"])
    
    # Determine start and end time for slicing
    t_min = 86400 * time_period["start"]
    t_max = 86400 * time_period["end"]

    # Find index associated with t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break

    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break

    # Reshape concentration matrices
    concentration_full = solution_full_model.y.reshape(num_tracers_full,num_boxes,len(solution_full_model.t))               # shape = (num_tracers, num_boxes, time)
    concentration_reduced = solution_reduced_model.y.reshape(num_tracers_reduced,num_boxes,len(solution_reduced_model.t))   # shape = (num_tracers, num_boxes, time)

    # Find target indices
    target_indices_full = []
    target_indices_reduced = []
    
    for name in tracer_names:
        # Index in solution_full_model
        if name in target_names: target_indices_full.append(tracer_names.index(name))
    
        # Split tracer name and constituent
        tracer,constituent = name.split("_")
    
        # Index in solution_reduced_model
        if tracer in tracer_map_reduced:
            if name in target_names:
                # Index of target constituent in target tracer
                const_index = tracers_reduced[tracer].composition.index(constituent)
                target_indices_reduced.append(tracer_map_reduced[tracer][const_index])

    # Extract concentration over the time period
    target_quality_full = np.zeros((len(target_indices_full),num_boxes,(index_t_max_full-index_t_min_full)))                # shape = (num_targets, num_boxes, time)
    target_quality_reduced = np.zeros((len(target_indices_reduced),num_boxes,(index_t_max_reduced-index_t_min_reduced)))    # shape = (num_targets, num_boxes, time)

    for index_full in range(index_t_min_full,index_t_max_full):
        target_quality_full[...,(index_full - index_t_min_full)] = concentration_full[target_indices_full,:,index_full]
    for index_reduced in range(index_t_min_reduced,index_t_max_reduced):
        target_quality_reduced[...,(index_reduced - index_t_min_reduced)] = concentration_reduced[target_indices_reduced,:,index_reduced]

    # Collapse target quality matrix
    collapsed_full = np.sum(target_quality_full, axis=0, keepdims=True)         # shape = (1, num_boxes, time)
    collapsed_reduced = np.sum(target_quality_reduced, axis=0, keepdims=True)   # shape = (1, num_boxes, time)

    # Find peak concentration
    peak_full = np.max(collapsed_full)
    peak_reduced = np.max(collapsed_reduced)

    # Compute the error with respect to full solution
    error = np.abs(100 * np.abs(peak_full - peak_reduced) / peak_full)

    return error, solution_reduced_model


def time_of_peak(t_span, time_period, base_element, physical, tracer_names, target_names, solution_full_model, c0_reduced, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names):
    """ calculates error in the time of peak concentration of target tracer (sum of targets if multiple are provided)
    """
    num_tracers_full = len(tracer_names)
    num_tracers_reduced = num_tracers(tracer_map_reduced)
    
    num_boxes = c0_reduced.shape[1]     # number of boxes in water column
    
    # solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names, False), 
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bgc_rate_eqns(time, base_element, conc, num_tracers_reduced, physical, tracers_reduced, tracer_map_reduced, tracer_type_reduced, indices_to_retain, removed_tracer_names), 
                                 t_span, c0_reduced.ravel(), method='RK23')#, max_step=physical["dt"])
    
    # Determine start and end time for slicing
    t_min = 86400 * time_period["start"]
    t_max = 86400 * time_period["end"]
    
    # Find index associated with t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
    
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # Reshape concentration matrices
    concentration_full = solution_full_model.y.reshape(num_tracers_full,num_boxes,len(solution_full_model.t))               # shape = (num_tracers, num_boxes, time)
    concentration_reduced = solution_reduced_model.y.reshape(num_tracers_reduced,num_boxes,len(solution_reduced_model.t))   # shape = (num_tracers, num_boxes, time)
    
    # Find target indices
    target_indices_full = []
    target_indices_reduced = []
    
    for name in tracer_names:
        # Index in solution_full_model
        if name in target_names: target_indices_full.append(tracer_names.index(name))
    
        # Split tracer name and constituent
        tracer,constituent = name.split("_")
    
        # Index in solution_reduced_model
        if tracer in tracer_map_reduced:
            if name in target_names:
                # Index of target constituent in target tracer
                const_index = tracers_reduced[tracer].composition.index(constituent)
                target_indices_reduced.append(tracer_map_reduced[tracer][const_index])
    
    # Extract concentration over the time period
    target_quality_full = np.zeros((len(target_indices_full),num_boxes,(index_t_max_full-index_t_min_full)))                # shape = (num_targets, num_boxes, time)
    target_quality_reduced = np.zeros((len(target_indices_reduced),num_boxes,(index_t_max_reduced-index_t_min_reduced)))    # shape = (num_targets, num_boxes, time)
    
    for index_full in range(index_t_min_full,index_t_max_full):
        target_quality_full[...,(index_full - index_t_min_full)] = concentration_full[target_indices_full,:,index_full]
    for index_reduced in range(index_t_min_reduced,index_t_max_reduced):
        target_quality_reduced[...,(index_reduced - index_t_min_reduced)] = concentration_reduced[target_indices_reduced,:,index_reduced]
    
    # Collapse target quality matrix
    collapsed_full = np.sum(target_quality_full, axis=0, keepdims=True)         # shape = (1, num_boxes, time)
    collapsed_reduced = np.sum(target_quality_reduced, axis=0, keepdims=True)   # shape = (1, num_boxes, time)

    # Slice list of time
    time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    time_reduced = solution_reduced_model[index_t_min_reduced:index_t_max_reduced]
    
    # Find the index of peak concentration
    peak_index_full = np.argmax(collapsed_full,axis=2)
    peak_index_reduced = np.argmax(collapsed_reduced,axis=2)

    # Find the time of peak concentration
    peak_time_full = time_full[peak_index_full]
    peak_time_reduced = time_reduced[peak_index_reduced]
    
    # Compute the error with respect to full solution
    error = np.zeros(num_boxes)
    for i in range(0,num_boxes):
        error[i] = np.abs(100 * np.abs(peak_time_full[i] - peak_time_reduced[i]) / peak_time_full[i])

    error = np.max(error)

    return error, solution_reduced_model
