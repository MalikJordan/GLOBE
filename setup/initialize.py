import json
import numpy as np
import os

from setup.functional_groups import Inorganic, LivingOrganic, NonLivingOrganic

class Tracer:
    """
    Description: Class for describing tracer properties.
    """
    def __init__(self, long_name, index):
        self.name = long_name   # long name
        self.conc = 0.          # concentration
        self.d_dt = 0.          # rate of change
        self.index = index      # index in concentration matrix


def coordinate_system(parameters):

    parameters["dz"] = parameters["column_depth"] / parameters["num_layers"]
    parameters["z"] = np.linspace(parameters["dz"]/2, parameters["column_depth"] - parameters["dz"]/2, parameters["num_layers"])

    return parameters


def import_model():
    
    # Open json file containing tracer information
    file_path = os.getcwd() + "/npzd-kidston.json"
    # file_path = os.getcwd() + "/model_description.json"
    with open(file_path) as read_model:
        model_info = json.load(read_model)
        model = model_info["tracers"]
        parameters = model_info["parameters"]
        reactions = model_info["reactions"]

    # Define coordinate system
    parameters["water_column"] = coordinate_system(parameters["water_column"])

    # Initialize dictionary of reactions for each tracer
    tracers = {}
    index = 0

    # Initialize array of concentrations
    conc = []

    # ----------------------------------------------------------------------------------------------------
    # Read Model Tracers
    # ----------------------------------------------------------------------------------------------------
    # Bacterioplankton
    if "bacterioplankton" in model:
        bacterioplankton = model["bacterioplankton"]
        bac = []
        for key in bacterioplankton:
            bac.append(key)
            tracers[key] = Tracer(bacterioplankton[key]["long_name"], index)
            conc.append(0.)
            index += 1

    else:
        print("Bacterioplankton not included in model.")
        bac = None

    # ----------------------------------------------------------------------------------------------------
    # Detritus
    if "detritus" in model:
        detritus = model["detritus"]
        det = []
        for key in detritus:
            det.append(key)
            tracers[key] = Tracer(detritus[key]["long_name"], index)
            conc.append(0.)
            index += 1

    else:
        print("Detritus not included in model.")
        det = None

    # ----------------------------------------------------------------------------------------------------
    # Inorganic Nutients
    if "inorganic" in model:
        inorganic = model["inorganic"]
        inorg = []
        for key in inorganic:
            inorg.append(key)
            tracers[key] = Tracer(inorganic[key]["long_name"], index)
            conc.append(0.)
            index += 1

    else:
        print("Inorganic nutrients not included in model.")
        inorg = None

    # ----------------------------------------------------------------------------------------------------
    # Phytoplankton
    if "phytoplankton" in model:
        phytoplankton = model["phytoplankton"]
        phyto = []
        for key in phytoplankton:
            phyto.append(key)
            tracers[key] = Tracer(phytoplankton[key]["long_name"], index)
            conc.append(0.)
            index += 1

    else:
        print("Phytoplankton not included in model.")
        phyto = None

    # ----------------------------------------------------------------------------------------------------
    # Zooplankton
    if "zooplankton" in model:
        zooplankton = model["zooplankton"]
        zoo = []
        for key in zooplankton:
            zoo.append(key)
            tracers[key] = Tracer(zooplankton[key]["long_name"], index)
            conc.append(0.)
            index += 1

    else:
        print("Zooplankton not included in model.")
        zoo = None

    return conc, parameters, tracers
