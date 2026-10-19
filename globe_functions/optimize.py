import copy
import os
import yaml
from pathlib import Path
# from optimization.build import extract_parameters
from optimization.shared import extract_parameters, reduce_parameter_space
from optimization.morris import run_morris_screening
from optimization.sobol import run_sobol_analysis
from optimization.optimize import run_optimization


def optimize_bgc_model(base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical):

    # if physical["simulation"]["reduce"] == True:    model_file_path = os.getcwd() + '/reduced_model.yaml' 
    # else:   model_file_path = os.getcwd() + '/model.yaml'
    if physical["simulation"]["reduce"] == True:    model_file_path = Path.cwd() / 'reduced_model.yaml' 
    else:   model_file_path = Path.cwd() / 'model.yaml'
    
    # Read bgc input file
    with open(model_file_path, 'r') as f:
        model_info = yaml.full_load(f)
        model = model_info["tracers"]
        reactions = model_info["reactions"]

    # Begin optimization
    if physical["simulation"]["optimize"] == True:
        # # Read dakota input file
        # dakota_file_path = os.getcwd() + '/dakota.yaml'
        # with open(dakota_file_path, 'r') as f:
        #     dakota_configuration = yaml.full_load(f)

        # optimization_file_path = os.getcwd() + '/optimization.yaml'
        optimization_file_path = Path.cwd() / 'optimization.yaml'
        with open(optimization_file_path, 'r') as f:
            optimization_configuration = yaml.full_load(f)

        parameter_bounds = optimization_configuration["parameter_bounds"]
        objective_targets = optimization_configuration["objective_targets"]
        morris_info = optimization_configuration["morris_screening"]
        sobol_info = optimization_configuration["sobol_SA"]
        opt_info = optimization_configuration["parameter_optimization"]
        
        # Extract parameters and calculate bounds
        # parameters, numeric_parameters, bounds = extract_parameters(base_element, model, dakota_configuration)
        parameters, independent_parameters, independent_bounds, constraints = extract_parameters(base_element, model, parameter_bounds)

        # Run Morris Screening on independent parameters
        num_trajectories = morris_info["trajectories"]
        num_levels = morris_info["levels"]
        num_workers = morris_info["workers"]
        threshold = morris_info["threshold"]
        # X, Y, results, retained, removed, sensitivity_scores = run_morris_screening(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_trajectories, num_levels, threshold, seed=1234)
        X, Y, results, retained, removed, sensitivity_scores = run_morris_screening(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_trajectories, num_levels, num_workers, threshold, seed=1234)

        # Reduce parameter space
        independent_parameters, independent_bounds, constraints = reduce_parameter_space(independent_parameters, independent_bounds, removed, retained, constraints)

        # Run Sobol' Sensitivity Analysis on remaining parameters
        num_samples = sobol_info["samples"]
        threshold = sobol_info["threshold"]
        num_workers = sobol_info["workers"]
        confidence_factor = sobol_info["confidence_factor"]
        # X, Y, results, retained, removed, relative_importance = run_sobol_analysis(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_samples, threshold, confidence_factor, seed=1234)
        X, Y, results, retained, removed, relative_importance = run_sobol_analysis(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, num_samples, threshold, confidence_factor, num_workers, seed=1234)

        # Reduce parameter space
        independent_parameters, independent_bounds, constraints = reduce_parameter_space(independent_parameters, independent_bounds, removed, retained, constraints)
        
        # Run parameter optimization
        max_iter = opt_info["max_iter"]
        tolerance = opt_info["tolerance"]
        num_workers = opt_info["workers"]
        optimized_parameters, result = run_optimization(base_element, model, reactions, physical, independent_parameters, independent_bounds, constraints, objective_targets, max_iter, tolerance, num_workers)
        x=1

    return base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical