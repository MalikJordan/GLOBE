import copy
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from scipy.io import loadmat
from globe_functions.simulation import run_simulation
from setup.initialize import import_bgc_model
from setup.other_functions import create_function_inputs
import pickle

_worker_context = None

def initialize_worker(base_element, model, reactions, physical, independent_parameters, constraints, objective_targets):
    global _worker_context

    _worker_context = {"base_element": base_element, "model": model, "reactions": reactions, "physical": physical, "independent_parameters": independent_parameters, "constraints": constraints, "objective_targets": objective_targets}

# ----------------------------------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------------------------------
def evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets):
    
    # Reconstruct all parameters, including dependent OM partitions
    sampled_parameters = reconstruct_parameters(independent_parameters, constraints, sample)

    # Update the tracer configuration
    # sampled_model = update_parameter_values(model, sampled_parameters)
    sampled_model = update_parameter_values(copy.deepcopy(model), sampled_parameters)

    # Initialize model with the sampled parameters
    tracers = import_bgc_model(base_element, sampled_model, reactions, physical)
    concentration, sinking, tracer_map, tracer_type = create_function_inputs(physical["simulation"]["iters"],tracers)

    # Run BGC simulation
    conc_daily, conc_monthly, npp_daily, npp_monthly = run_simulation(base_element, concentration, sinking, tracers, tracer_map, tracer_type, physical)

    # Calculate objective function
    objective_target_data = get_objective_target_data(conc_monthly[...,:-1], npp_monthly[...,:-1], objective_targets, tracers, tracer_map)
    objective = calculate_objective_function(objective_target_data)

    return objective


# def evaluate_one_sample(args):
#     """
#     Evaluate a single parameter sample
#     """
#     # Extract arguments
#     (base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets,) = args

#     # Evaluate sample and return the objective function
#     return evaluate_sample(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets,)


def evaluate_one_sample(sample):
    """
    Evaluate a single parameter sample
    """
    # Extract arguments
    ctx = _worker_context

    # Evaluate sample and return the objective function
    return evaluate_sample(ctx["base_element"], ctx["model"], ctx["reactions"], ctx["physical"], ctx["independent_parameters"], ctx["constraints"], sample, ctx["objective_targets"],)


# def parallel_evaluation(X, base_element, model, reactions, physical, independent_parameters, constraints, objective_targets, num_workers=4):
#     """
#     Evaluate samples in parallel
#     """
#     objects = {
#         "model": model,
#         "reactions": reactions,
#         "physical": physical,
#         "independent_parameters": independent_parameters,
#         "constraints": constraints,
#     }

#     for name, obj in objects.items():
#         try:
#             pickle.dumps(obj)
#             print(f"{name}: pickleable")
#         except Exception as exc:
#             print(f"{name}: NOT pickleable — {exc}")

#     # Create taks for the given samples
#     tasks = [(base_element, model, reactions, physical, independent_parameters, constraints, sample, objective_targets,) for sample in X]

#     # Run parallel evaluation
#     with ProcessPoolExecutor(max_workers=num_workers) as executor:    Y = list(executor.map(evaluate_one_sample, tasks))

#     return np.asarray(Y, dtype=float)


def parallel_evaluation(X, base_element, model, reactions, physical, independent_parameters, constraints, objective_targets, num_workers=4):
    """
    Evaluate samples in parallel
    """

    # Run parallel evaluation
    with ProcessPoolExecutor(max_workers=num_workers, initializer=initialize_worker, initargs=(base_element, model, reactions, physical, independent_parameters, constraints, objective_targets)) as executor:    
        # Initialize configuration for each worker to reduce overhead
        Y = list(executor.map(evaluate_one_sample, X))

    return np.asarray(Y, dtype=float)



# ----------------------------------------------------------------------------------------------------
# Parameter Functions
# ----------------------------------------------------------------------------------------------------
def extract_parameters(base_element, model, parameter_bounds):
    """
    Definition:: extracts all numeric parameters for every tracer and calculates their bounds for Morris Screening, Sobol' Sensitivity Analysis, and parameter optimization.
    
    ----------------------------------------------------------------------
    Paramters
    ----------
    base_element : str
        bale element used by the model (e.g. "c", "n", "p")
    model : dict
        Model configuration containing tracer parameter dictionaries
    parameter_bounds : dict
        Dictionary containing lower and upper bound multipliers

    ----------------------------------------------------------------------
    Returns
    ----------
    parameters : dict
        All flattened tracer parameters, including non-numeric values
    independent_parameters : dict
        Numeric parameters excluding dependent OM partition parameters
    independent_bounds : dict
        Bounds for independent parameters, including feasible bounds for OM partition parameters
    constraints : list
        OM partition constraints
    """
    # Extract and flatten parameters for all tracers
    parameters = {}
    for tracer_name, tracer in model.items():
        tracer_parameters = tracer.get("parameters",{})
        parameters.update(flatten_parameters(tracer_parameters,tracer_name))

    # Retain all parameter types but only optimize numeric values
    numeric_parameters = {name: value for name,value in parameters.items() if isinstance(value,(int,float)) and not isinstance(value,bool)}

    # Remove "base_temp", cell quotas for base element, pH, and OCMIP data from numeric parameter dictionary
    for key in list(numeric_parameters.keys()):
        if "base_temp" in key:  del numeric_parameters[key]
        elif "pH" in key:       del numeric_parameters[key]
        elif "OCMIP" in key:    del numeric_parameters[key]
        elif "cell_quota" in key and key.endswith("." + base_element):  del numeric_parameters[key]

    # Set constraints for om_partitiions
    constraints = set_om_partition_constraints(numeric_parameters)

    # Read multipliers for parameter bounds
    lower = parameter_bounds["lower"]
    upper = parameter_bounds["upper"]

    # Verify realistic bounds
    if not (0.0 <= lower <= 1.0 <= upper):  raise ValueError("Parameter-bound multipliers must satisfy 0 <= lower_bound <= 1 <= upper_bound.")
    
    # Calculate bounds using constant values from Dakota YAML file
    bounds = {}
    for name,value in numeric_parameters.items():
        # Use multipliers if initial_value > 0. else apply explorable range
        if value > 0.:  lower_value , upper_value = lower*value , upper*value
        else:   lower_value , upper_value = 0. , 0.1

        # Bound om_partition [0,1]
        if ".om_partition." in name: lower_value, upper_value = 0., 1.

        bounds[name] = {"initial": value, "lower": lower_value, "upper": upper_value}

    # Remove one dependent parameter from each om_partition group.
    dependent_parameters = {constraint["parameters"][-1] for constraint in constraints}

    # Set independent parameters andd bounds for future perturbations (Morris screening, Sobol' sensitivity analysis, and Dakota parameter optimmization)
    independent_parameters = {name: value for name, value in numeric_parameters.items() if name not in dependent_parameters}
    independent_bounds = {name: bound for name, bound in bounds.items() if name not in dependent_parameters}

    # Calculate feasible individual bounds for partition variables.
    # These bounds are necessary, but do not alone enforce the
    # joint sum constraint.
    for constraint in constraints:
        names = constraint["parameters"]
        total = constraint["total"]

        for name in names[:-1]:
            other_names = [other for other in names if other != name]

            # feasible_lower = max(bounds[name]["lower"],total - sum(bounds[other]["upper"] for other in other_names),)
            feasible_lower = max(bounds[name]["lower"], total - sum(bounds[other]["upper"] for other in other_names))
            # feasible_upper = min(bounds[name]["upper"], total - sum(bounds[other]["lower"] for other in other_names),)
            feasible_upper = min(bounds[name]["upper"], total - sum(bounds[other]["lower"] for other in other_names))

            if feasible_lower > feasible_upper: raise ValueError(f"Infeasible OM partition bounds for {name}")

            independent_bounds[name]["lower"] = feasible_lower
            independent_bounds[name]["upper"] = feasible_upper

    return parameters, independent_parameters, independent_bounds, constraints


    # return parameters, numeric_parameters, bounds


def reconstruct_parameters(independent_parameters, constraints, sample):
    """
    Definition:: Reconstruct the complete parameter dictionary, including dependent OM partition parameters, for one sensitivity-analysis sample.

    ----------------------------------------------------------------------
    Parameters
    ----------
    independent_parameters : dict
        Independent parameters and their initial values
    constraints : list
        OM partition constraints returned by extract_parameters()
    sample : array-like
        Sampled values in the order of independent_parameters.keys()

    ----------------------------------------------------------------------
    Returns
    -------
    parameters : dict
        Complete sampled parameter dictionary, including dependent values
    """
    # Create list of parameter names
    parameter_names = list(independent_parameters.keys())

    # Verify correct number of parameters
    if len(parameter_names) != len(sample): raise ValueError("The number of independent parameters must match the sample length.")

    # Assign sampled values to independent parameters.
    parameters = {name: float(value) for name,value in zip(parameter_names, sample)}

    # Reconstruct each dependent OM partition parameter.
    for constraint in constraints:
        names = constraint["parameters"]
        total = constraint["total"]
        dependent_name = names[-1]

        dependent_value = total - sum(parameters[name] for name in names[:-1])

        # Value for partition cannot be zero
        if dependent_value < 0.:    raise ValueError(f"Invalid OM partition sample: dependent parameter '{dependent_name}' would be {dependent_value}.")

        # Avoid tiny negative values from floating-point roundoff.
        parameters[dependent_name] = max(0.0, dependent_value)

    return parameters


def reduce_parameter_space(independent_parameters, independent_bounds, removed, retained, constraints):

    # Parse constraints list for removed parameters
    constraints = [constraint for constraint in constraints if any(param in retained for param in constraint["parameters"])]

    # Reduce independent parameters and bounds dictionaries
    for key in removed:
        if key in list(independent_parameters.keys()):  
            del independent_parameters[key]
            del independent_bounds[key]

    return independent_parameters, independent_bounds, constraints


def update_parameter_values(tracers, parameters):
    """
    Definition:: Update nested YAML parameters using Dakota parameter names

    Example:
        'phyto.growth.max_photo_rate' --> tracers['phyto']['parameters']['growth']['max_photo_rate']
    """
    for parameter_name,value in parameters.items():
        # Split parameters into parts
        parts = parameter_name.split(".")

        # Get traer name and path to parameter
        tracer_name = parts[0]
        parameter_path = parts[1:]

        # Navigate to tracer's parameter dictionary
        target = tracers[tracer_name].setdefault("parameters",{})

        # Navigate through nested parameter dictionaries
        for key in parameter_path[:-1]: target = target.setdefault(key,{})

        # Update final parameter value
        target[parameter_path[-1]] = float(value)

    return tracers

# ----------------------------------------------------------------------------------------------------
# Parameter Constraints for om_partition
# ----------------------------------------------------------------------------------------------------
def set_om_partition_constraints(numeric_parameters):
    
    groups = defaultdict(list)

    for key in numeric_parameters:
        if ".om_partition." not in key: continue

        # Extract source and destination from parameter key
        source, destination = key.split(".om_partition.",1)

        # Identify element from destination (destination = tracer.element)
        element = destination.rsplit(".",1)[-1]
        groups[(source,element)].append(key)

    # Set constraints to ensure sum of destinations = 1 for each source
    constraints = []
    for (source,element),names in groups.items():
        if len(names) > 1:
            constraints.append({"group": source, "element": element, "parameters": names, "total": 1.})

    return constraints


def apply_om_partition_constraints(model, sampled, constraints):
    """
    Definition:: applies constraints to organic matter partition to ensure all destination groups for an individual source group sum to 1
    """

    applied = dict(sampled)
    constraints = constraints or []

    # Apply independent values first
    for name,value in applied.items():  set_dotted_value(model, name, value)

    # Reconstruct dependent fractions
    for constraint in constraints:
        names = constraint.get("parameters",[])
        if len(names) < 2:  continue
        total = float(constraint.get("total",1.))
        dependent_name = names[-1]
        independent_names = names[:-1]
        dependent_value = total - sum( float(get_dotted_value(model,name)) for name in independent_names)

        if dependent_value < 0:
            raise ValueError(f"Invalid sum constraint {constraint!r}: dependent fraction would be {dependent_value:.8g}")
        dependent_value = max(0.,dependent_value)
        set_dotted_value(model, dependent_name, dependent_value)
        applied[dependent_name] = dependent_value

    return applied

# ----------------------------------------------------------------------------------------------------
# Parameter Mapping
# ----------------------------------------------------------------------------------------------------
def flatten_parameters(data,prefix=""):
    """
    Definition:: flattens nested dictionaries into dotted parameter paths
    """
    parameters = {}
    for name,value in data.items():
        path = f"{prefix}.{name}" if prefix else name

        if isinstance(value,dict):  parameters.update(flatten_parameters(value,path))
        else:   parameters[path] = value

    return parameters


def get_dotted_value(map, dotted_path):
    current: any = map
    for part in dotted_path.split("."): current = current[part]

    return current


def set_dotted_value(map, dotted_path, value):
    """
    Definition:: set a value using a dotted dictionary path
    """
    parts = dotted_path.split(".")
    current: any = map

    for part in parts[:-1]:
        if not isinstance(current, dict) or part not in current:
            raise KeyError(f"Cannot resolve parameter path {dotted_path!r} at {part!r}")
        current = current[part]
    if not isinstance(current, dict) or parts[-1] not in current:
        raise KeyError(f"Cannot resolve parameter path {dotted_path!r}")
    current[parts[-1]] = float(value)

# ----------------------------------------------------------------------------------------------------
# Objective Functions
# ----------------------------------------------------------------------------------------------------
def calculate_objective_function(objective_target_data):

    # Extract targets from target data
    objective_targets = list(objective_target_data.keys())

    # Read observational data
    observations = get_observational_data()
    
    # Calculate normalization scales for observational data
    scales = calulate_objective_scales(observations, objective_targets)

    nrmse = []
    target_rmse = {}

    # Trim observational data for rmse calculations
    for target in objective_targets:
        # Extract data and trim observational data
        model_values = objective_target_data[target]
        num_months = model_values.shape[1]
        observed_values = observations[target][:,:num_months]

        # Verify validity of data
        valid = (np.isfinite(model_values) & np.isfinite(observed_values))
        if not np.any(valid):   raise ValueError(f"{target}: no matching valid monthly values.")

        # Calculate rmse
        rmse = np.sqrt( np.mean(model_values[valid] - observed_values[valid]) ** 2 )
        target_rmse[target] = float(rmse)

        # Append nrmse
        nrmse.append(rmse/scales[target])

    # Calculate objective (average nrmse of targets)
    objective = float(np.mean(nrmse))

    return objective


def calulate_objective_scales(observations, objective_targets):
    """
    Definition:: Calculates fixed normalization scales from monthly observations.

    Parameters
    ----------
    observations : dict
        Monthly observational averages, keyed by target name.
    objective_targets : list
        Targets included in the objective function.

    Returns
    -------
    scales : dict
        Sample standard deviation for each target.
    """

    # Calculate standard deviation for each target
    scales = {}
    for target in objective_targets:
        values = np.asarray(observations[target], dtype=np.float64)
        values = values[np.isfinite(values)]

        # if values.size < 2:
        #     raise ValueError(
        #         f"{target} requires at least two valid monthly observations."
        #     )

        scale = np.std(values, ddof=1)
        if scale <= 0:  raise ValueError(f"{target} has zero observational standard deviation.")

        scales[target] = float(scale)

    return scales


def get_observational_data():

    # Read observational data
    observations = {}
    # path = os.getcwd() + '/optimization/Data_ML/'
    path = Path.cwd() / 'optimization/Data_ML/'
    
    # Load Data
    chl_data = loadmat(path / 'Chla_1yr_climatology.mat')
    oxy_data  = loadmat(path / 'Oxy_1yr_climatology.mat')
    nit_data  = loadmat(path / 'Nitrate_1yr_climatology.mat')
    phos_data = loadmat(path / 'Phos_1yr_climatology.mat')
    pon_data  = loadmat(path / 'PON_1yr_climatology.mat')
    npp_data  = loadmat(path / 'NPP_1yr_climatology.mat')
    
    # Extract Upper 150m
    observations["chl"] = chl_data['Var'][:150,:]
    observations["o2"]  = oxy_data['Var'][:150,:]
    observations["no3"]  = nit_data['Var'][:150,:]
    observations["po4"] = phos_data['Var'][:150,:]
    observations["pon"]  = pon_data['Var'][:150,:]
    observations["npp"]  = npp_data['Var'][:150,:]

    return observations


def get_objective_target_data(concentration, net_primary_production, objective_targets, tracers, tracer_map):
    """
    Definition:: Extracts simulated data for objective function targets
    
    Parameters
    ----------
    concentration : array
        Simulated monthly averages for all tracers
    net_primary_production : array
        Simulated monthly average for net primary production
    objective_targets : list[str]
        Targets included in the objective function
    tracers : dict
        Tracers present in model
    tracer_map : dict
        Map of tracer indices in concentration array. ke = tracer name , val= list of indices
    
    Returns
    -------
    objective_target_data : dict
        Simulated monthly averages for all objective targets
    """

    objective_target_data = {}

    if "o2" in objective_targets:   objective_target_data["o2"] = concentration[tracer_map["o2"][0]]
    if "no3" in objective_targets:  objective_target_data["no3"] = concentration[tracer_map["no3"][0]]
    if "po4" in objective_targets:  objective_target_data["po4"] = concentration[tracer_map["po4"][0]]
    if "npp" in objective_targets:  objective_target_data["npp"] = net_primary_production
    if "chl" in objective_targets:
        # Initialize array for chlorophyll data
        chlorophyll = np.zeros_like(concentration[0,...])
        for trac in list(tracers.keys()):
            # Locate phytoplankton groups
            if tracers[trac].type == "phytoplankton":
                # Find chlorophyll index
                if 'chl' in tracers[trac].composition:  chl_index = tracers[trac].composition.index('chl')
                else:   continue

                # Add chlorophyll concentration of phytoplankton group to data array
                chlorophyll += concentration[tracer_map[trac][chl_index]]

        objective_target_data["chl"] = chlorophyll

    if "pon" in objective_targets:
        # Initialize array for chlorophyll data
        pon = np.zeros_like(concentration[0,...])
        for trac in list(tracers.keys()):
            # Locate phytoplankton groups
            if tracers[trac].type in ["phytoplankton","zooplankton"]:
                # Find nitrogen index
                if 'n' in tracers[trac].composition:  n_index = tracers[trac].composition.index('n')
                else:   continue
        
                # Add nitrogen concentration of phytoplankton group to data array
                pon += concentration[tracer_map[trac][n_index]]

            elif tracers[trac].type == "detritus":
                if tracers[trac].form == "particulate":
                    # Find nitrogen index
                    if 'n' in tracers[trac].composition:  n_index = tracers[trac].composition.index('n')
                    else:   continue

                    # Add nitrogen concentration of phytoplankton group to data array
                    pon += concentration[tracer_map[trac][n_index]]
        
        objective_target_data["pon"] = pon
        
    return objective_target_data

