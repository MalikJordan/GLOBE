# import dakota.interfacing as di

import numpy as np
import os
import yaml
from collections import defaultdict
from pathlib import Path


# def build_environment(env_config):
#     lines = ["environment"]

#     if "tabular_data_file" in env_config:
#         lines.append("  tabular_data")  # 2 spaces before "tabular_data" (don't change)
#         lines.append(f"    tabular_data_file = '{env_config['tabular_data_file']}'")  # 4 spaces before "tabular_data_file" (don't change)

#     return "\n".join(lines) + "\n\n"


# def build_interface(intf_config):
#     lines = ["interface"]

#     if intf_config.get("fork", True):
#         lines.append("  fork")
#     else:
#         lines.append("  system")

#     for key,value in intf_config.items():
#         if key == "fork":
#             continue

#         if isinstance(value, bool):
#             if value:
#                 lines.append(f"    {key}")
#         else:
#             lines.append(f"    {key} = '{value}'" if key else f"    {key} = {value}")

#     return "\n".join(lines) + "\n\n"


# def build_method(mthd_config):
#     method_type = mthd_config.get("type", "sampling")

#     lines = ["method", f"  {method_type}"]  # 2 spaces before "{method_type}" (don't change)

#     for key,value in mthd_config.items():
#         if key == "type":
#             continue

#         if isinstance(value, bool):
#             if value:
#                 lines.append(f"    {key}")  # 4 spaces before "{key}" (don't change)
#         else:
#             lines.append(f"    {key} = {value}")    # 4 spaces before "{key}" (don't change)
        
#     return "\n".join(lines) + "\n"


# # def build_method_integer_controls(mthd_config):
# #     method_type = mthd_config.get("type", "sampling")

# #     DI_MTHD = METHOD_MAP[method_type]

# #     # Method-specific integer inputs
# #     if method_type == "sampling":
# #         DI_MI = [
# #             mthd_config.get("samples", 100),    # default to 100
# #             mthd_config.get("seed", 12345)      # default to 12345
# #         ]
# #     elif method_type in ["moga", "nsga2"]:
# #         DI_MI = [
# #             mthd_config.get("population_size", 100),    # default to 100
# #             mthd_config.get("max_generations", 100)     # default to 100
# #         ]
# #     elif method_type == "nlpql_sqp":
# #         DI_MI = [
# #             mthd_config.get("max_iterations", 100)      # default to 100
# #         ]
# #     elif method_type == "bayes_calibration":
# #         DI_MI = [
# #             mthd_config.get("chain_samples", 1000)      # default to 1000
# #         ]
# #     else:
# #         raise ValueError(f"Unsupported method: {method_type}")
    
# #     return DI_MTHD, DI_MI


# def build_model(mod_config):
#     model_type = mod_config.get("type", "single")

#     # 2 spaces before "{model_type}" (don't change)
#     return f"""
# model
#   {model_type}\n
# """


# def build_response(resp_config):
#     lines = ["responses"]

#     n_obj = resp_config.get("num_objective_functions", 1)   # default to 1 objective function
#     lines.append(f"  num_objective_functions = {n_obj}")    # 2 spaces before "num_objective_functions" (don't change)

#     gradients = resp_config.get("gradients", False) # default to no_gradients
#     hessians = resp_config.get("hessians", False)   # default to no_hessians

#     if gradients:
#         lines.append("  analytic_gradients")    # 2 spaces before "analytic_gradients" (don't change)
#     else:
#         lines.append("  no_gradients")  # 2 spaces before "no_gradients" (don't change)

#     if hessians:
#         lines.append("  analytic_hessians") # 2 spaces before "analytic_hessians" (don't change)
#     else:
#         lines.append("  no_hessians")   # 2 spaces before "no_hessians" (don't change)

#     return "\n".join(lines) + "\n\n"


# def build_variables(params, vars_config=None):  # keep vars_config to potentially extended model to include discrete variables for optimization
#     n = len(params)

#     # 1 space (don't change)
#     lower_bounds = " ".join(str(p["lower"]) for p in params)
#     upper_bounds = " ".join(str(p["upper"]) for p in params)
#     initial_point = " ".join(str(p["initial"]) for p in params)
#     descriptors = " ".join(f"'{p['name']}'" for p in params)

#     lines = [
#         "variables",    # no spaces (don't change)
#         f"  continuous_design = {n}",   # 2 spaces before "variables" (don't change)
#         f"    initial_point     {initial_point}", # 4 spaces before "initial_point" (don't change)
#         f"    lower_bounds      {lower_bounds}",  # 4 spaces before "lower_bounds" (don't change)
#         f"    upper_bounds      {upper_bounds}",  # 4 spaces before "upper_bounds" (don't change)
#         f"    descriptors       {descriptors}"    # 4 spaces before "descriptors" (don't change)
#     ]
    
#     return "\n".join(lines) + "\n\n"

# # def extract_parameters(parameter_set, prefix=""):

# #     params = [] # Create empty list to store parameters

# #     for key,value in parameter_set.items(): # Parse through parameter_set dictionary
# #         new_key = f"{prefix}.{key}" if prefix else key  # Define new key in case it contains additional levels

# #         # Check if this is a parameter leaf
# #         if isinstance(value, dict) and all(k in value for k in ["lower", "initial", "upper"]):
# #         # Only contains "lower", "initial", and "upper" == parameter --> add to list of parameters
# #             params.append({
# #                 "name": new_key,
# #                 "lower": value["lower"],
# #                 "initial": value["initial"],
# #                 "upper": value["upper"]
# #             })
# #         elif isinstance(value, dict):
# #         # Contains additional 
# #             params.extend(extract_parameters(value, new_key))

# #     return params

def extract_parameters(base_element, model, dakota_configuration):
    """
    Definition:: extracts all numeric parameters for every tracer
    """
    parameters = {}
    for tracer_name, tracer in model.items():
        tracer_parameters = tracer.get("parameters",{})
        parameters.update(flatten_parameters(tracer_parameters,tracer_name))

    # Retain all parameter types but only optimize numeric values
    numeric_parameters = {name: value for name,value in parameters.items() if isinstance(value,(int,float)) and not isinstance(value,bool)}

    # Remove "base_temp" and cell quotas for base element from numeric parameter dictionary
    for key in list(numeric_parameters.keys()):
        if "base_temp" in key:  del numeric_parameters[key]
        elif "cell_quota" in key and key.endswith("." + base_element):  del numeric_parameters[key]

    # Set constraints for om_partitiions
    constraints = set_om_partition_constraints(numeric_parameters)

    # Calculate bounds using constant values from Dakota YAML file
    lower = dakota_configuration["parameter_bounds"]["lower"]
    upper = dakota_configuration["parameter_bounds"]["upper"]
    bounds = {}

    for name,value in numeric_parameters.items():
        lower_value = value * lower
        upper_value = value * upper

        if lower_value > upper_value:   lower_value,upper_value = upper_value,lower_value
        elif upper_value == 0.:         upper_value = 0.1
        # if lower_value == upper_value:  raise ValueError(f"Parameter '{name}' has identical bounds. A nonzero range is required.")

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


def update_yaml_parameters(tracers, parameters):
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


def read_dakota_parameters(parameters_file, results_file):
    """
    Read Dakota parameter values and prepare its response object.

    Returns
    -------
    parameters : dict
        Parameter identifiers mapped to their trial values.
    dakota_results : Dakota Results object
        Object used to write the QoI back to Dakota.
    """
    params, dakota_results = di.read_parameters_file(
        parameters_file=str(parameters_file),
        results_file=str(results_file),
    )

    parameters = {
        name: float(params[name])
        for name in params
    }

    return parameters, dakota_results


def dakota_driver(parameters_file, results_file, reduced=True):

    # Select model configuration
    yaml_file = Path.cwd() / ('reduced_model.yaml' if reduced else 'model.yaml')

    # Load base configuration
    with yaml_file.open('r') as file:
        model_info = yaml.safe_load(file)
        base_element = model_info["base_element"]
        tracers = model_info["tracers"]
        reactions = model_info["reactions"]

    # Read parameter values from dakota file
    # parameters = read_dakota_parameters(parameters_file)

    # Update YAML parameter values
    # update_yaml_parameters(tracers,parameters)

    pass