import yaml
import os
import sys
from setup.initialize import coordinate_system
from setup.bacteria import Bacteria
from setup.detritus import Detritus
from setup.inorganic import Inorganic
from setup.phytoplankton import Phytoplankton
from setup.zooplankton import MesoZooplankton, MicroZooplankton

file = os.getcwd() + '/npzd.yaml'
with open(file, 'r') as f:
    model_info = yaml.full_load(f)
    model = model_info["tracers"]
    parameters = model_info["parameters"]
    reactions = model_info["reactions"]

# Update coordinate system
parameters["water_column"] = coordinate_system(parameters["water_column"])

tracers = {}
index = 0

for key in model:
    if model[key]["type"] == "bacteria":
        tracers[key] = Bacteria(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    elif model[key]["type"] == "detritus":
        tracers[key] = Detritus(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    elif model[key]["type"] == "inorganic":
        tracers[key] = Inorganic(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    elif model[key]["type"] == "phytoplankton":
        tracers[key] = Phytoplankton(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    elif model[key]["type"] == "mesozooplankton":
        tracers[key] = MesoZooplankton(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    elif model[key]["type"] == "microzooplankton":
        tracers[key] = MicroZooplankton(index, model[key]["long_name"], model[key]["composition"], model[key]["ic"])
    else:
        sys.exit("Warning: Functional group '" + model[key]["type"] + "' not accepted. Please review documentation and make necessary changes.")
        
# Add necessary components to Phytoplankton and Zooplankton groups
for key in reactions:
    if key["type"] == "grazing":    # Add prey to zooplankton (used in rate calculations to determine sum of grazing rates)
        tracers[key["produced"]].add_prey(key["consumed"])
    if key["type"] == "uptake":     # Add nutrient to phytoplankton (used in rate calculations for nutrient limitation)
        tracers[key["produced"]].add_nutrient(key["consumed"])

print()


# def import_model():
    
#     # Open json file containing tracer information
#     file_path = os.getcwd() + "/npzd-riley.json"
#     # file_path = os.getcwd() + "/model_description.json"
#     with open(file_path) as read_model:
#         model_info = json.load(read_model)
#         model = model_info["tracers"]
#         parameters = model_info["parameters"]
#         reactions = model_info["reactions"]

#     # Define coordinate system
#     parameters["water_column"] = coordinate_system(parameters["water_column"])

#     # Initialize dictionary of reactions for each tracer
#     tracers = {}
#     index = 0

#     # Initialize array of concentrations
#     conc = []

#     # ----------------------------------------------------------------------------------------------------
#     # Read Model Tracers
#     # ----------------------------------------------------------------------------------------------------
#     # Bacterioplankton
#     if "bacterioplankton" in model:
#         bacterioplankton = model["bacterioplankton"]
#         bac = []
#         for key in bacterioplankton:
#             bac.append(key)
#             tracers[key] = Tracer(bacterioplankton[key]["long_name"], index)
#             conc.append(0.)
#             index += 1

#     else:
#         print("Bacterioplankton not included in model.")
#         bac = None

#     # ----------------------------------------------------------------------------------------------------
#     # Detritus
#     if "detritus" in model:
#         detritus = model["detritus"]
#         det = []
#         for key in detritus:
#             det.append(key)
#             tracers[key] = Tracer(detritus[key]["long_name"], index)
#             conc.append(0.)
#             index += 1

#     else:
#         print("Detritus not included in model.")
#         det = None

#     # ----------------------------------------------------------------------------------------------------
#     # Inorganic Nutients
#     if "inorganic" in model:
#         inorganic = model["inorganic"]
#         inorg = []
#         for key in inorganic:
#             inorg.append(key)
#             tracers[key] = Tracer(inorganic[key]["long_name"], index)
#             conc.append(0.)
#             index += 1

#     else:
#         print("Inorganic nutrients not included in model.")
#         inorg = None

#     # ----------------------------------------------------------------------------------------------------
#     # Phytoplankton
#     if "phytoplankton" in model:
#         phytoplankton = model["phytoplankton"]
#         phyto = []
#         for key in phytoplankton:
#             phyto.append(key)
#             tracers[key] = Tracer(phytoplankton[key]["long_name"], index)
#             conc.append(0.)
#             index += 1

#     else:
#         print("Phytoplankton not included in model.")
#         phyto = None

#     # ----------------------------------------------------------------------------------------------------
#     # Zooplankton
#     if "zooplankton" in model:
#         zooplankton = model["zooplankton"]
#         zoo = []
#         for key in zooplankton:
#             zoo.append(key)
#             tracers[key] = Tracer(zooplankton[key]["long_name"], index)
#             conc.append(0.)
#             index += 1

#     else:
#         print("Zooplankton not included in model.")
#         zoo = None

#     return conc, parameters, tracers
# #