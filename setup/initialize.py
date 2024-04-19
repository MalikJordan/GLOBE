import json
import os
from setup.functional_groups import Bacterioplankton, Inorganic, OrganicMatter, Phytoplankton, Zooplankton

def initialize_model_tracers():
    """
    Description: Initializes BGC model. Creates a dictionary of all tracers present in model
                 and classes for each group. 
    """

    # ----------------------------------------------------------------------------------------------------
    num_boxes = 1
    # Initialize dictionary of tracers present in model
    tracers = {}
    count_tracers = 0

    # Open json file containing tracer information
    file_path = os.getcwd() + "/tracers.json"
    with open(file_path) as read_tracers:
        tracer_info = json.load(read_tracers)

    # ----------------------------------------------------------------------------------------------------
    # Add non-living tracers (inorganic, dissolved organic matter, particulate organic matter)

    if "inorganic" in tracer_info:
        inorganic = tracer_info["inorganic"]
        inorg = Inorganic(inorganic, num_boxes)
        for key in inorganic:
            if "constituents" in inorganic:
                constituents = inorganic[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1            
    else:
        print("Inrganic tracers not listed in " + file_path)
        inorg = {}

    if "dissolved_organic_matter" in tracer_info:
        dissolved_inorganic_matter = tracer_info["dissolved_organic_matter"]
        dom = {}
        for key in dissolved_inorganic_matter:
            dom[key] = OrganicMatter(dissolved_inorganic_matter[key]["long_name"], dissolved_inorganic_matter[key]["constituents"], num_boxes)
            if "constituents" in dissolved_inorganic_matter[key]:
                constituents = dissolved_inorganic_matter[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Dissolved Organic Matter not listed in " + file_path)
        dom = {}

    if "particulate_organic_matter" in tracer_info:
        particulate_organic_matter = tracer_info["particulate_organic_matter"]
        pom = {}
        for key in particulate_organic_matter:
            pom[key] = OrganicMatter(particulate_organic_matter[key]["long_name"], particulate_organic_matter[key]["constituents"], num_boxes)
            if "constituents" in particulate_organic_matter[key]:
                constituents = particulate_organic_matter[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Particulate Organic Matter not listed in " + file_path)
        pom = {}

    # ----------------------------------------------------------------------------------------------------
    # Add living tracers (bacterioplankton, phytoplankton, zooplankton)

    if "bacterioplankton" in tracer_info:
        bacterioplankton = tracer_info["bacterioplankton"]
        bac = {}
        for key in bacterioplankton:
            bac[key] = Bacterioplankton(bacterioplankton[key]["long_name"], num_boxes)
            if "constituents" in bacterioplankton[key]:
                constituents = bacterioplankton[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Bacterioplankton not listed in " + file_path)
        bac = {}

    if "phytoplankton" in tracer_info:
        phytoplankton = tracer_info["phytoplankton"]
        phyto = {}
        for key in phytoplankton:
            phyto[key] = Phytoplankton(phytoplankton[key]["long_name"], phytoplankton[key]["constituents"], num_boxes)
            if "constituents" in phytoplankton[key]:
                constituents = phytoplankton[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Phytoplankton not listed in " + file_path)
        phyto = {}

    if "zooplankton" in tracer_info:
        zooplankton = tracer_info["zooplankton"]
        zoo = {}
        for key in zooplankton:
            zoo[key] = Zooplankton(zooplankton[key]["long_name"], num_boxes)
            if "constituents" in zooplankton[key]:
                constituents = zooplankton[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Zooplankton not listed in " + file_path)
        zoo = {}

    return tracers, bac, dom, inorg, phyto, pom, zoo