import numpy as np
import os
import yaml


def build_environment(env_config):
    lines = ["environment"]

    if "tabular_data_file" in env_config:
        lines.append("  tabular_data")  # 2 spaces before "tabular_data" (don't change)
        lines.append(f"    tabular_data_file = '{env_config['tabular_data_file']}'")  # 4 spaces before "tabular_data_file" (don't change)

    return "\n".join(lines) + "\n\n"


def build_interface(intf_config):
    lines = ["interface"]

    if intf_config.get("fork", True):
        lines.append("  fork")
    else:
        lines.append("  system")

    for key,value in intf_config.items():
        if key == "fork":
            continue

        if isinstance(value, bool):
            if value:
                lines.append(f"    {key}")
        else:
            lines.append(f"    {key} = '{value}'" if key else f"    {key} = {value}")

    return "\n".join(lines) + "\n\n"


def build_method(mthd_config):
    method_type = mthd_config.get("type", "sampling")

    lines = ["method", f"  {method_type}"]  # 2 spaces before "{method_type}" (don't change)

    for key,value in mthd_config.items():
        if key == "type":
            continue

        if isinstance(value, bool):
            if value:
                lines.append(f"    {key}")  # 4 spaces before "{key}" (don't change)
        else:
            lines.append(f"    {key} = {value}")    # 4 spaces before "{key}" (don't change)
        
    return "\n".join(lines) + "\n"


# def build_method_integer_controls(mthd_config):
#     method_type = mthd_config.get("type", "sampling")

#     DI_MTHD = METHOD_MAP[method_type]

#     # Method-specific integer inputs
#     if method_type == "sampling":
#         DI_MI = [
#             mthd_config.get("samples", 100),    # default to 100
#             mthd_config.get("seed", 12345)      # default to 12345
#         ]
#     elif method_type in ["moga", "nsga2"]:
#         DI_MI = [
#             mthd_config.get("population_size", 100),    # default to 100
#             mthd_config.get("max_generations", 100)     # default to 100
#         ]
#     elif method_type == "nlpql_sqp":
#         DI_MI = [
#             mthd_config.get("max_iterations", 100)      # default to 100
#         ]
#     elif method_type == "bayes_calibration":
#         DI_MI = [
#             mthd_config.get("chain_samples", 1000)      # default to 1000
#         ]
#     else:
#         raise ValueError(f"Unsupported method: {method_type}")
    
#     return DI_MTHD, DI_MI


def build_model(mod_config):
    model_type = mod_config.get("type", "single")

    # 2 spaces before "{model_type}" (don't change)
    return f"""
model
  {model_type}\n
"""


def build_response(resp_config):
    lines = ["responses"]

    n_obj = resp_config.get("num_objective_functions", 1)   # default to 1 objective function
    lines.append(f"  num_objective_functions = {n_obj}")    # 2 spaces before "num_objective_functions" (don't change)

    gradients = resp_config.get("gradients", False) # default to no_gradients
    hessians = resp_config.get("hessians", False)   # default to no_hessians

    if gradients:
        lines.append("  analytic_gradients")    # 2 spaces before "analytic_gradients" (don't change)
    else:
        lines.append("  no_gradients")  # 2 spaces before "no_gradients" (don't change)

    if hessians:
        lines.append("  analytic_hessians") # 2 spaces before "analytic_hessians" (don't change)
    else:
        lines.append("  no_hessians")   # 2 spaces before "no_hessians" (don't change)

    return "\n".join(lines) + "\n\n"


def build_variables(params, vars_config=None):  # keep vars_config to potentially extended model to include discrete variables for optimization
    n = len(params)

    # 1 space (don't change)
    lower_bounds = " ".join(str(p["lower"]) for p in params)
    upper_bounds = " ".join(str(p["upper"]) for p in params)
    initial_point = " ".join(str(p["initial"]) for p in params)
    descriptors = " ".join(f"'{p['name']}'" for p in params)

    lines = [
        "variables",    # no spaces (don't change)
        f"  continuous_design = {n}",   # 2 spaces before "variables" (don't change)
        f"    initial_point     {initial_point}", # 4 spaces before "initial_point" (don't change)
        f"    lower_bounds      {lower_bounds}",  # 4 spaces before "lower_bounds" (don't change)
        f"    upper_bounds      {upper_bounds}",  # 4 spaces before "upper_bounds" (don't change)
        f"    descriptors       {descriptors}"    # 4 spaces before "descriptors" (don't change)
    ]
    
    return "\n".join(lines) + "\n\n"

def extract_parameters(parameter_set, prefix=""):

    params = [] # Create empty list to store parameters

    for key,value in parameter_set.items(): # Parse through parameter_set dictionary
        new_key = f"{prefix}.{key}" if prefix else key  # Define new key in case it contains additional levels

        # Check if this is a parameter leaf
        if isinstance(value, dict) and all(k in value for k in ["lower", "initial", "upper"]):
        # Only contains "lower", "initial", and "upper" == parameter --> add to list of parameters
            params.append({
                "name": new_key,
                "lower": value["lower"],
                "initial": value["initial"],
                "upper": value["upper"]
            })
        elif isinstance(value, dict):
        # Contains additional 
            params.extend(extract_parameters(value, new_key))

    return params


def generate_dakota_input(yaml_file, output_file="dakota.in"):

    with open(yaml_file, "r") as f:
        config = yaml.safe_load(f)

    dakota_config = config.get("dakota", {})
    params_config = config.get("parameters", {})

    param_list = extract_parameters(params_config)

    text = ""
    text += build_environment(dakota_config.get("environment", {}))
    text += build_method(dakota_config.get("method", {}))
    text += build_model(dakota_config.get("model", {}))
    text += build_variables(param_list, dakota_config.get("variables", {}))
    text += build_interface(dakota_config.get("interface", {}))
    text += build_response(dakota_config.get("responses", {}))
    
    with open(output_file, "w") as f:
        f.write(text)

    return param_list


# METHOD_MAP = {
#     "sampling": 1,
#     "moga": 2,
#     "ngsa2": 3,
#     "nlpql_sqp": 4,
#     "bayes_calibration": 5
# }
# with open("dakota.yaml", "r") as f:
#     config = yaml.safe_load(f)

# param_list = extract_parameters(config["parameters"])

# num_params = len(param_list)
# print("Number of continuous variables:", num_params)

# dakota_vars_block = build_variables(param_list)

# with open("dakota.in", "w") as f:
#     f.write(dakota_vars_block)

generate_dakota_input("dakota.yaml")
x=1