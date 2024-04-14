import json
import os


def initialize_model_tracers():

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

    if "dissolved_organic_matter" in tracer_info:
        dom = tracer_info["dissolved_organic_matter"]
        for key in dom:
            if "constituents" in dom[key]:
                constituents = dom[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Dissolved Organic Matter not listed in " + file_path)

    if "particulate_organic_matter" in tracer_info:
        pom = tracer_info["particulate_organic_matter"]
        for key in pom:
            if "constituents" in pom[key]:
                constituents = pom[key]["constituents"]
                for const in constituents:
                    tracers[count_tracers] = key + '_' + const
                    count_tracers += 1
            else:
                tracers[count_tracers] = key
                count_tracers += 1
    else:
        print("Particulate Organic Matter not listed in " + file_path)

    # ----------------------------------------------------------------------------------------------------
    # Add living tracers (bacterioplankton, phytoplankton, zooplankton)

    if "bacterioplankton" in tracer_info:
        bacterioplankton = tracer_info["bacterioplankton"]
        for key in bacterioplankton:
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

    if "phytoplankton" in tracer_info:
        phytoplankton = tracer_info["phytoplankton"]
        for key in phytoplankton:
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

    if "zooplankton" in tracer_info:
        zooplankton = tracer_info["zooplankton"]
        for key in zooplankton:
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

    return tracers