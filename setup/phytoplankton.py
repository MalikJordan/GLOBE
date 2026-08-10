import os
import sys
import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutrient_limitation, monod, temperature_dependence, tracer_elements, switch
from fractions import Fraction
from pom.check_phy import nitr_lim, phos_lim, mult_lim, photo
np.set_printoptions(precision=20)
class Phytoplankton():
    """
    
    """

    # def __init__(self, abbrev, base_element, iters, num_layers, reactions, **tracer):
    def __init__(self, abbrev, base_element, physical, reactions, **tracer):
        
        # Variales that will be used later ---------------------------------------------------------------
        num_layers = physical["water_column"]["num_layers"]
        iters = physical["simulation"]["iters"]
        composition = physical["initial_conditions"][abbrev]["composition"]     # Initial concentrations
        if "scale" in physical["initial_conditions"][abbrev]:   scale = physical["initial_conditions"][abbrev]["scale"]     # Scaling factor (if initial concentration is split between multiple phytoplankton groups)
        else:   scale = 1.  # No scaling

        # Add important keys ---------------------------------------------------------------
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Nutrient limitation
        self.nutrient_limitation = tracer["parameters"]["nutrient_limitation"]
        self.nutrient_limitation_factor = Dict.empty(key_type=types.unicode_type,value_type=types.ListType(float64[:]))
        self.nutrient_colimitation_factor = 0.
        for key in self.nutrient_limitation:
            # add "nh4_inhibited" to nitrate limitation parameters if not in dictionary (needed for uptake)
            if key == "no3" and "nh4_inhibited" not in self.nutrient_limitation["no3"]:   self.nutrient_limitation["no3"]["nh4_inhibited"] = False
            if key == "nh4" and "nh4_inhibited" not in self.nutrient_limitation["nh4"]:   self.nutrient_limitation["nh4"]["nh4_inhibited"] = False
        
        # Exudation
        if "exudation" in tracer["parameters"]:
            # Create float option numbers for use in numba typed.List
            if "method" in tracer["parameters"]["exudation"]:
                if tracer["parameters"]["exudation"]["method"] == "constant":           tracer["parameters"]["exudation"]["method"] = 1
                elif tracer["parameters"]["exudation"]["method"] == "photosynthesis":   tracer["parameters"]["exudation"]["method"] = 2
                elif tracer["parameters"]["exudation"]["method"] == "uptake":           tracer["parameters"]["exudation"]["method"] = 3

            self.exudation_ids = []
            self.exudation_params = []
            
            for key,val in tracer["parameters"]["exudation"].items():
                self.exudation_ids.append(key)
                # if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.exudation_params.append(np.float64(val))

        # Growth
        if "growth" in tracer["parameters"]:
            # Create float option numbers for use in numba typed.List
            if "light_limitation" in tracer["parameters"]["growth"]:
                if tracer["parameters"]["growth"]["light_limitation"] == "monod":       tracer["parameters"]["growth"]["light_limitation"] = -1
                elif tracer["parameters"]["growth"]["light_limitation"] == "geider":    tracer["parameters"]["growth"]["light_limitation"] = -2
                elif tracer["parameters"]["growth"]["light_limitation"] == "platt":     tracer["parameters"]["growth"]["light_limitation"] = -3
                elif tracer["parameters"]["growth"]["light_limitation"] == "smith":     tracer["parameters"]["growth"]["light_limitation"] = -4

            if "light_location" in tracer["parameters"]["growth"]:
                if tracer["parameters"]["growth"]["light_location"] == "top":           tracer["parameters"]["growth"]["light_location"] = 1
                elif tracer["parameters"]["growth"]["light_location"] == "middle":      tracer["parameters"]["growth"]["light_location"] = 2
                elif tracer["parameters"]["growth"]["light_location"] == "integrated":  tracer["parameters"]["growth"]["light_location"] = 3

            if "silicate_limitation" in tracer["parameters"]["growth"]:
                if tracer["parameters"]["growth"]["silicate_limitation"] == True:       tracer["parameters"]["growth"]["silicate_limitation"] = 1
                else:                                                                   tracer["parameters"]["growth"]["silicate_limitation"] = 0
            else:   tracer["parameters"]["growth"]["silicate_limitation"] = 0
            
            # Define numeric option key for variable max growth rate
            if isinstance(tracer["parameters"]["growth"]["max_photo_rate"],str):        tracer["parameters"]["growth"]["max_photo_rate"] = -1
            
            # Define numeric option key for type growth rate type (if max photo rate is calculated using eppley formation)
            if "type" in tracer["parameters"]["growth"]:
                if tracer["parameters"]["growth"]["type"] == "base_b":      tracer["parameters"]["growth"]["type"] = 1
                elif tracer["parameters"]["growth"]["type"] == "standard":  tracer["parameters"]["growth"]["type"] = 2

            # Translate from fraction string to float64 (if necessary)
            if "convert_o2" in tracer["parameters"]["growth"] and isinstance(tracer["parameters"]["growth"]["convert_o2"],str):
                tracer["parameters"]["growth"]["convert_o2"] = np.float64(Fraction(tracer["parameters"]["growth"]["convert_o2"]))

            self.growth_ids = List.empty_list(unicode_type)
            self.growth_params = List.empty_list(float64)

            for key,val in tracer["parameters"]["growth"].items():
                self.growth_ids.append(key)
                # if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.growth_params.append(np.float64(val))

        # Lysis
        if "lysis" in tracer["parameters"]:
            # Create typed.Dict for dissolved/particulate apportioning factor
            self.lysis_apportioning_factor = Dict.empty(
                key_type=types.unicode_type,
                value_type=types.ListType(types.unicode_type)
            )

            if "apportioning_factor" in tracer["parameters"]["lysis"]:
                for key,val in tracer["parameters"]["lysis"]["apportioning_factor"].items():    
                    val_list = List.empty_list(unicode_type)
                    for const in val:   val_list.append(const)
                    self.lysis_apportioning_factor[key] = val_list

            # Create float option numbers for use in numba typed.List
            if "method" in tracer["parameters"]["lysis"]:
                if tracer["parameters"]["lysis"]["method"] == "cell_quota":     tracer["parameters"]["lysis"]["method"] = 1
                elif tracer["parameters"]["lysis"]["method"] == "constant":     tracer["parameters"]["lysis"]["method"] = 2

                if "convert_lysis" in tracer["parameters"]["lysis"] and isinstance(tracer["parameters"]["lysis"]["convert_lysis"],str):
                    tracer["parameters"]["lysis"]["convert_lysis"] = np.array([Fraction(tracer["parameters"]["lysis"]["convert_lysis"])],dtype=np.float64)

            self.lysis_ids = List.empty_list(unicode_type)
            self.lysis_params = List.empty_list(float64[:])

            for key,val in tracer["parameters"]["lysis"].items():
                if key != "apportioning_factor":
                    self.lysis_ids.append(key)
                    if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    self.lysis_params.append(val)

        # Respiration
        if "respiration" in tracer["parameters"]:
            # Translate from fraction string to float64 (if necessary)
            if "convert_o2" in tracer["parameters"]["respiration"]:
                if isinstance(tracer["parameters"]["respiration"]["convert_o2"],str):
                    tracer["parameters"]["respiration"]["convert_o2"] = np.array([Fraction(tracer["parameters"]["respiration"]["convert_o2"])],dtype=np.float64)
                else:
                    tracer["parameters"]["respiration"]["convert_o2"] = np.array([tracer["parameters"]["respiration"]["convert_o2"]],dtype=np.float64)
            if "convert_co2" in tracer["parameters"]["respiration"]:
                if isinstance(tracer["parameters"]["respiration"]["convert_co2"],str):
                    tracer["parameters"]["respiration"]["convert_co2"] = np.array([Fraction(tracer["parameters"]["respiration"]["convert_co2"])],dtype=np.float64)
                else:
                    tracer["parameters"]["respiration"]["convert_co2"] = np.array([tracer["parameters"]["respiration"]["convert_co2"]],dtype=np.float64)
            
            self.respiration_ids = List.empty_list(unicode_type)
            self.respiration_params = List.empty_list(float64[:])

            for key,val in tracer["parameters"]["respiration"].items():
                self.respiration_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.respiration_params.append(val)

        # Uptake
        if "uptake" in tracer["parameters"]:
            # Create typed.Dict of uptake parameters and ids
            self.uptake_ids = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(unicode_type))
            self.uptake_params = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(float64))

            # Potential coupled uptake "link" ids
            # [1] no3, [2] nh4, [3] po4, [4] fe, [5] sio4
            coupled_uptake_codes = {"no3": np.float64(1.), "nh4": np.float64(2.), "po4": np.float64(3.), "fe": np.float64(4.), "sio4": np.float64(5.)}
            
            # Nested typed.Dict for coupled uptake
            self.coupled_uptake = Dict.empty(
                key_type=types.unicode_type, 
                value_type=types.DictType(types.unicode_type, types.ListType(types.float64))
            )

            # Nitrate
            if "no3" in tracer["parameters"]["uptake"]:
                # Create float option numbers for use in numba typed.List
                if "basis" in tracer["parameters"]["uptake"]["no3"]:
                    if tracer["parameters"]["uptake"]["no3"]["basis"] == "constant":      tracer["parameters"]["uptake"]["no3"]["basis"] = 1
                    elif tracer["parameters"]["uptake"]["no3"]["basis"] == "growth":      tracer["parameters"]["uptake"]["no3"]["basis"] = 2
                    elif tracer["parameters"]["uptake"]["no3"]["basis"] == "nutrient":    tracer["parameters"]["uptake"]["no3"]["basis"] = 3

                if "strategy" in tracer["parameters"]["uptake"]["no3"]:
                    if tracer["parameters"]["uptake"]["no3"]["strategy"] == "independent":    tracer["parameters"]["uptake"]["no3"]["strategy"] = 1
                    elif tracer["parameters"]["uptake"]["no3"]["strategy"] == "coupled":
                        tracer["parameters"]["uptake"]["no3"]["strategy"] = 2

                        # Create inner typed.Dict for coupled uptake
                        coupled_uptake_no3 = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.float64))

                        # links = linked nutrient(s) to use uptake rates of
                        # Have to change list items to option codes before applying to numba typed.Dict
                        links_list = List.empty_list(float64)   # Initialize list of coupled uptake option ids
                        for key in coupled_uptake_codes:    # Append list
                            if key in tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["link"]:  links_list.append(coupled_uptake_codes[key])
                        coupled_uptake_no3["links"] = links_list  # Add to coupled uptake dictionary
                        
                        # method = option code for method of linked calculation
                        method_list = List.empty_list(float64)
                        if len(links_list) > 1:
                            if tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["method"] == "max":    coupled_uptake_no3["method"] = 1
                            elif tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["method"] == "min":  coupled_uptake_no3["method"] = 2
                            elif tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["method"] == "sum":  coupled_uptake_no3["method"] = 3
                            elif tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["method"] == "product":  coupled_uptake_no3["method"] = 4
                        else:
                            coupled_uptake_no3["method"] = 1

                        # Translate from fraction string to float64 (if necessary)
                        if isinstance(tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["convert_uptake"],str):
                            frac = Fraction(tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["convert_uptake"])
                            frac = np.float64(frac)
                            coupled_uptake_no3["convert_uptake"] = List(frac)
                        else:
                            coupled_uptake_no3["convert_uptake"] = np.float64(tracer["parameters"]["uptake"]["no3"]["coupled_uptake"]["convert_uptake"])

                        # Add to overall coupled uptake typed.Dict
                        self.coupled_uptake["no3"] = coupled_uptake_no3

                if "form" in tracer["parameters"]["uptake"]["no3"]:
                    if tracer["parameters"]["uptake"]["no3"]["form"] == "affinity":       tracer["parameters"]["uptake"]["no3"]["form"] = 1
                    elif tracer["parameters"]["uptake"]["no3"]["form"] == "constituent":  tracer["parameters"]["uptake"]["no3"]["form"] = 2

                # Add uptake keys,values to numba typed.Lists
                uptake_no3_ids = List.empty_list(unicode_type)
                uptake_no3_params = List.empty_list(float64)
                for key,val in tracer["parameters"]["uptake"]["no3"].items():
                    if key == "coupled_uptake": continue
                    uptake_no3_ids.append(key)
                    if not isinstance(val, np.float64): val = np.float64(val)  # Convert type to array of floats for typed.List
                    uptake_no3_params.append(val)

                # Add uptake_n lists to numba typed.Dicts
                self.uptake_ids["no3"] = uptake_no3_ids
                self.uptake_params["no3"] = uptake_no3_params

            # Ammonium
            if "nh4" in tracer["parameters"]["uptake"]:
                # Create float option numbers for use in numba typed.List
                if "basis" in tracer["parameters"]["uptake"]["nh4"]:
                    if tracer["parameters"]["uptake"]["nh4"]["basis"] == "constant":      tracer["parameters"]["uptake"]["nh4"]["basis"] = 1
                    elif tracer["parameters"]["uptake"]["nh4"]["basis"] == "growth":      tracer["parameters"]["uptake"]["nh4"]["basis"] = 2
                    elif tracer["parameters"]["uptake"]["nh4"]["basis"] == "nutrient":    tracer["parameters"]["uptake"]["nh4"]["basis"] = 3

                if "strategy" in tracer["parameters"]["uptake"]["nh4"]:
                    if tracer["parameters"]["uptake"]["nh4"]["strategy"] == "independent":    tracer["parameters"]["uptake"]["nh4"]["strategy"] = 1
                    elif tracer["parameters"]["uptake"]["nh4"]["strategy"] == "coupled":
                        tracer["parameters"]["uptake"]["nh4"]["strategy"] = 2

                        # Create inner typed.Dict for coupled uptake
                        coupled_uptake_nh4 = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.float64))

                        # links = linked nutrient(s) to use uptake rates of
                        # Have to change list items to option codes before applying to numba typed.Dict
                        links_list = List.empty_list(float64)   # Initialize list of coupled uptake option ids
                        for key in coupled_uptake_codes:    # Append list
                            if key in tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["link"]:  links_list.append(coupled_uptake_codes[key])
                        coupled_uptake_nh4["links"] = links_list  # Add to coupled uptake dictionary
                        
                        # method = option code for method of linked calculation
                        method_list = List.empty_list(float64)
                        if len(links_list) > 1:
                            if tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["method"] == "max":    coupled_uptake_nh4["method"] = 1
                            elif tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["method"] == "min":  coupled_uptake_nh4["method"] = 2
                            elif tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["method"] == "sum":  coupled_uptake_nh4["method"] = 3
                            elif tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["method"] == "product":  coupled_uptake_nh4["method"] = 4
                        else:
                            coupled_uptake_nh4["method"] = 1

                        # Translate from fraction string to float64 (if necessary)
                        if isinstance(tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["convert_uptake"],str):
                            frac = Fraction(tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["convert_uptake"])
                            frac = np.float64(frac)
                            coupled_uptake_nh4["convert_uptake"] = List(frac)
                        else:
                            coupled_uptake_nh4["convert_uptake"] = np.float64(tracer["parameters"]["uptake"]["nh4"]["coupled_uptake"]["convert_uptake"])

                        # Add to overall coupled uptake typed.Dict
                        self.coupled_uptake["nh4"] = coupled_uptake_nh4

                if "form" in tracer["parameters"]["uptake"]["nh4"]:
                    if tracer["parameters"]["uptake"]["nh4"]["form"] == "affinity":       tracer["parameters"]["uptake"]["nh4"]["form"] = 1
                    elif tracer["parameters"]["uptake"]["nh4"]["form"] == "constituent":  tracer["parameters"]["uptake"]["nh4"]["form"] = 2

                # Add uptake keys,values to numba typed.Lists
                uptake_nh4_ids = List.empty_list(unicode_type)
                uptake_nh4_params = List.empty_list(float64)
                for key,val in tracer["parameters"]["uptake"]["nh4"].items():
                    if key == "coupled_uptake": continue
                    uptake_nh4_ids.append(key)
                    if not isinstance(val, np.float64): val = np.float64(val)  # Convert type to array of floats for typed.List
                    uptake_nh4_params.append(val)

                # Add uptake_n lists to numba typed.Dicts
                self.uptake_ids["nh4"] = uptake_nh4_ids
                self.uptake_params["nh4"] = uptake_nh4_params

            # Phosphate
            if "po4" in tracer["parameters"]["uptake"]:
                # Create float option numbers for use in numba typed.List
                if "basis" in tracer["parameters"]["uptake"]["po4"]:
                    if tracer["parameters"]["uptake"]["po4"]["basis"] == "constant":      tracer["parameters"]["uptake"]["po4"]["basis"] = 1
                    elif tracer["parameters"]["uptake"]["po4"]["basis"] == "growth":      tracer["parameters"]["uptake"]["po4"]["basis"] = 2
                    elif tracer["parameters"]["uptake"]["po4"]["basis"] == "nutrient":    tracer["parameters"]["uptake"]["po4"]["basis"] = 3

                if "strategy" in tracer["parameters"]["uptake"]["po4"]:
                    if tracer["parameters"]["uptake"]["po4"]["strategy"] == "independent":    tracer["parameters"]["uptake"]["po4"]["strategy"] = 1
                    elif tracer["parameters"]["uptake"]["po4"]["strategy"] == "coupled":
                        tracer["parameters"]["uptake"]["po4"]["strategy"] = 2

                        # Create inner typed.Dict for coupled uptake
                        coupled_uptake_po4 = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.float64))

                        # links = linked nutrient(s) to use uptake rates of
                        # Have to change list items to option codes before applying to numba typed.Dict
                        links_list = List.empty_list(float64)   # Initialize list of coupled uptake option ids
                        for key in coupled_uptake_codes:    # Append list
                            if key in tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["link"]:  links_list.append(coupled_uptake_codes[key])
                        coupled_uptake_po4["links"] = links_list  # Add to coupled uptake dictionary
                        
                        # method = option code for method of linked calculation
                        method_list = List.empty_list(float64)
                        if len(links_list) > 1:
                            if tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["method"] == "max":        coupled_uptake_po4["method"] = 1
                            elif tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["method"] == "min":      coupled_uptake_po4["method"] = 2
                            elif tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["method"] == "sum":      coupled_uptake_po4["method"] = 3
                            elif tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["method"] == "product":  coupled_uptake_po4["method"] = 4
                        else:
                            method_list.append(1.)  # type doesn't matter if uptake is only linked to one nutrient

                        coupled_uptake_po4["method"] = method_list
                        # Translate from fraction string to float64 (if necessary)
                        if isinstance(tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["convert_uptake"],str):
                            frac = Fraction(tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["convert_uptake"])
                            frac = np.float64(frac)
                            coupled_uptake_po4["convert_uptake"] = List(frac)
                        else:
                            coupled_uptake_po4["convert_uptake"] = np.float64(tracer["parameters"]["uptake"]["po4"]["coupled_uptake"]["convert_uptake"])

                        # Add to overall coupled uptake typed.Dict
                        self.coupled_uptake["po4"] = coupled_uptake_po4

                if "form" in tracer["parameters"]["uptake"]["po4"]:
                    if tracer["parameters"]["uptake"]["po4"]["form"] == "affinity":       tracer["parameters"]["uptake"]["po4"]["form"] = 1
                    elif tracer["parameters"]["uptake"]["po4"]["form"] == "constituent":  tracer["parameters"]["uptake"]["po4"]["form"] = 2

                # Add uptake keys,values to numba typed.Lists
                uptake_po4_ids = List.empty_list(unicode_type)
                uptake_po4_params = List.empty_list(float64)
                for key,val in tracer["parameters"]["uptake"]["po4"].items():
                    if key == "coupled_uptake": continue
                    uptake_po4_ids.append(key)
                    # if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    if not isinstance(val, np.float64): val = np.float64(val)  # Convert type to array of floats for typed.List
                    uptake_po4_params.append(val)

                # Add uptake_p lists to numba typed.Dicts
                self.uptake_ids["po4"] = uptake_po4_ids
                self.uptake_params["po4"] = uptake_po4_params

            # Iron
            if "fe" in tracer["parameters"]["uptake"]:
                # Create float option numbers for use in numba typed.List
                if "basis" in tracer["parameters"]["uptake"]["fe"]:
                    if tracer["parameters"]["uptake"]["fe"]["basis"] == "constant":     tracer["parameters"]["uptake"]["fe"]["basis"] = 1
                    elif tracer["parameters"]["uptake"]["fe"]["basis"] == "growth":     tracer["parameters"]["uptake"]["fe"]["basis"] = 2
                    elif tracer["parameters"]["uptake"]["fe"]["basis"] == "nutrient":   tracer["parameters"]["uptake"]["fe"]["basis"] = 3

                if "strategy" in tracer["parameters"]["uptake"]["fe"]:
                    if tracer["parameters"]["uptake"]["fe"]["strategy"] == "independent":   tracer["parameters"]["uptake"]["fe"]["strategy"] = 1
                    elif tracer["parameters"]["uptake"]["fe"]["strategy"] == "coupled":
                        tracer["parameters"]["uptake"]["fe"]["strategy"] = 2

                        # Create inner typed.Dict for coupled uptake
                        coupled_uptake_fe = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.float64))

                        # links = linked nutrient(s) to use uptake rates of
                        # Have to change list items to option codes before applying to numba typed.Dict
                        links_list = List.empty_list(float64)   # Initialize list of coupled uptake option ids
                        for key in coupled_uptake_codes:    # Append list
                            if key in tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["link"]:  links_list.append(coupled_uptake_codes[key])
                        coupled_uptake_fe["links"] = links_list  # Add to coupled uptake dictionary
                        
                        # method = option code for method of linked calculation
                        method_list = List.empty_list(float64)
                        if tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["method"] == "max":       coupled_uptake_fe["method"] = 1
                        elif tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["method"] == "min":     coupled_uptake_fe["method"] = 2
                        elif tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["method"] == "sum":     coupled_uptake_fe["method"] = 3
                        elif tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["method"] == "product": coupled_uptake_fe["method"] = 4

                        # Translate from fraction string to float64 (if necessary)
                        if isinstance(tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["convert_uptake"],str):
                            frac = Fraction(tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["convert_uptake"])
                            frac = np.float64(frac)
                            coupled_uptake_fe["conversion"] = List(frac)
                        else:
                            coupled_uptake_fe["conversion"] = np.float64(tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["convert_uptake"])

                        # Add to overall coupled uptake typed.Dict
                        self.coupled_uptake["fe"] = coupled_uptake_fe

                if "form" in tracer["parameters"]["uptake"]["fe"]:
                    if tracer["parameters"]["uptake"]["fe"]["form"] == "affinity":      tracer["parameters"]["uptake"]["fe"]["form"] = 1
                    elif tracer["parameters"]["uptake"]["fe"]["form"] == "constituent": tracer["parameters"]["uptake"]["fe"]["form"] = 2

                # Add uptake keys,values to numba typed.Lists
                uptake_fe_ids = List.empty_list(unicode_type)
                uptake_fe_params = List.empty_list(float64)
                for key,val in tracer["parameters"]["uptake"]["fe"].items():
                    uptake_fe_ids.append(key)
                    # if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    if not isinstance(val, np.float64): val = np.float64(val)  # Convert type
                    uptake_fe_params.append(val)

                # Add uptake_fe lists to numba typed.Dicts
                self.uptake_ids["fe"] = uptake_fe_ids
                self.uptake_params["fe"] = uptake_fe_params

            # Silicate
            if "sio4" in tracer["parameters"]["uptake"]:
                # Create float option numbers for use in numba typed.List
                if "basis" in tracer["parameters"]["uptake"]["sio4"]:
                    if tracer["parameters"]["uptake"]["sio4"]["basis"] == "constant":     tracer["parameters"]["uptake"]["sio4"]["basis"] = 1
                    elif tracer["parameters"]["uptake"]["sio4"]["basis"] == "growth":     tracer["parameters"]["uptake"]["sio4"]["basis"] = 2
                    elif tracer["parameters"]["uptake"]["sio4"]["basis"] == "nutrient":   tracer["parameters"]["uptake"]["sio4"]["basis"] = 3

                if "strategy" in tracer["parameters"]["uptake"]["sio4"]:
                    if tracer["parameters"]["uptake"]["sio4"]["strategy"] == "independent":   tracer["parameters"]["uptake"]["sio4"]["strategy"] = 1
                    elif tracer["parameters"]["uptake"]["sio4"]["strategy"] == "coupled":
                        tracer["parameters"]["uptake"]["sio4"]["strategy"] = 2

                        # Create inner typed.Dict for coupled uptake
                        coupled_uptake_sio4 = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.float64))

                        # links_n = linked nutrient(s) to use uptake rates of
                        # Have to change list items to option codes before applying to numba typed.Dict
                        links_list = List.empty_list(float64)   # Initialize list of coupled uptake option ids
                        for key in coupled_uptake_codes:    # Append list
                            if key in tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["link"]:  links_list.append(coupled_uptake_codes[key])
                        coupled_uptake_sio4["links"] = links_list  # Add to coupled uptake dictionary
                        
                        # links_n_method = option code for method of linked calculation
                        method_list = List.empty_list(float64)
                        if tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["method"] == "max":       coupled_uptake_sio4["method"] = 1
                        elif tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["method"] == "min":     coupled_uptake_sio4["method"] = 2
                        elif tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["method"] == "sum":     coupled_uptake_sio4["method"] = 3
                        elif tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["method"] == "product": coupled_uptake_sio4["method"] = 4

                        # Translate from fraction string to float64 (if necessary)
                        if isinstance(tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["convert_uptake"],str):
                            frac = Fraction(tracer["parameters"]["uptake"]["fe"]["coupled_uptake"]["convert_uptake"])
                            frac = np.float64(frac)
                            coupled_uptake_sio4["conversion"] = np.float64(frac)
                        else:
                            coupled_uptake_sio4["conversion"] = np.float64(tracer["parameters"]["uptake"]["sio4"]["coupled_uptake"]["convert_uptake"])

                        # Add to overall coupled uptake typed.Dict
                        self.coupled_uptake["sio4"] = coupled_uptake_sio4

                if "form" in tracer["parameters"]["uptake"]["sio4"]:
                    if tracer["parameters"]["uptake"]["sio4"]["form"] == "affinity":      tracer["parameters"]["uptake"]["sio4"]["form"] = 1
                    elif tracer["parameters"]["uptake"]["sio4"]["form"] == "constituent": tracer["parameters"]["uptake"]["sio4"]["form"] = 2

                # Add uptake keys,values to numba typed.Lists
                uptake_sio4_ids = List.empty_list(unicode_type)
                uptake_sio4_params = List.empty_list(float64)
                for key,val in tracer["parameters"]["uptake"]["sio4"].items():
                    uptake_sio4_ids.append(key)
                    # if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    if not isinstance(val, np.float64): val = np.float64(val)  # Convert type
                    uptake_sio4_params.append(val)

                # Add uptake_n lists to numba typed.Dicts
                self.uptake_ids["sio4"] = uptake_sio4_ids
                self.uptake_params["sio4"] = uptake_sio4_params

        # Temperature regulation
        if "temperature_regulation" in tracer["parameters"]:    
            # Initialize temperature regulation factor to array of 1. (if phytoplankton is temperature limited this will be updated later otherwise will stay as 1.)
            if num_layers > 1:  self.temp_regulation_factor = np.ones(num_layers-1, dtype=np.float64)
            else:   self.temp_regulation_factor = np.float64(1.)

            if "temp_limited" in tracer["parameters"]["temperature_regulation"]:    self.temp_limited = tracer["parameters"]["temperature_regulation"]["temp_limited"]
            else:   self.temp_limited = False

            if self.temp_limited:
                # Create float option numbers for use in numba typed.List
                if "function" in tracer["parameters"]["temperature_regulation"]:
                    if tracer["parameters"]["temperature_regulation"]["function"] == "arrhenius":   tracer["parameters"]["temperature_regulation"]["function"] = 1
                    elif tracer["parameters"]["temperature_regulation"]["function"] == "eppley":    tracer["parameters"]["temperature_regulation"]["function"] = 2
                    elif tracer["parameters"]["temperature_regulation"]["function"] == "q10":       tracer["parameters"]["temperature_regulation"]["function"] = 3

                # Add temperature regulation factor to dictionary
                # tracer["parameters"]["temperature_regulation"]["temp_regulation_factor"] = np.empty((0,),dtype=np.float64)

                self.temp_reg_ids = List.empty_list(unicode_type)
                self.temp_reg_params = List.empty_list(float64[:])
                
                for key,val in tracer["parameters"]["temperature_regulation"].items():
                    self.temp_reg_ids.append(key)
                    if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    self.temp_reg_params.append(val)
        else: self.temp_limited = False

        # Sedimentation
        if "sedimentation" in tracer["parameters"]:
            # Create attribute for background sinking velocity
            if tracer["parameters"]["sedimentation"]["sinking"] == True:
                if num_layers > 1:
                    self.sinking_velocity = np.ones(num_layers-1,dtype=np.float64) * tracer["parameters"]["sedimentation"]["background_sinking_rate"]
                    self.sinking_velocity[-1] = np.float64(tracer["parameters"]["sedimentation"]["burial_velocity"])
                else:
                    self.sinking_velocity = np.array(tracer["parameters"]["sedimentation"]["background_sinking_rate"],dtype=np.float64)
            else:   
                if num_layers > 1:  self.sinking_velocity = np.zeros(num_layers-1,dtype=np.float64)
                else:               self.sinking_velocity = np.array([0.],dtype=np.float64)
            
            # Add attribute for background sinking rate
            if "background_sinking_rate" in tracer["parameters"]["sedimentation"]:  self.background_sinking_rate = tracer["parameters"]["sedimentation"]["background_sinking_rate"]
            else:   self.background_sinking_rate = 0.

            # Add attribute for maximum sinking rate
            if "max_sinking_rate" in tracer["parameters"]["sedimentation"]: self.max_sinking_rate = tracer["parameters"]["sedimentation"]["max_sinking_rate"]
            else:   self.max_sinking_rate = 0.

            # Add attribute for sinking threshold
            if "sinking_threshold" in tracer["parameters"]["sedimentation"]: self.sinking_threshold = tracer["parameters"]["sedimentation"]["sinking_threshold"]
            else:   self.sinking_threshold = 0.

        # Add concentrations ---------------------------------------------------------------
        self.composition = List.empty_list(unicode_type)
        conc = []
        if len(composition) < 1:
            sys.exit("Phytoplankton: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            # Reorder "composition" so "base_element" is at the start of the list
            for key in list(composition):
                if key != base_element: composition[key] = composition.pop(key)

            for key in list(composition):
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    # Add constituent to composition/concentration (if not already in list)
                    if key in self.composition: pass
                    else:
                        self.composition.append(key)

                        # Set initial conditions
                        if isinstance(composition[key], str): # Read initial conditions from file
                            conc.append( np.fromfile(os.getcwd() + composition[key]) )

                        elif isinstance(composition[key], (int,float)): # Create array of initial conditions based off initial value
                            if num_layers == 1: # 0d configuration
                                conc.append( composition[key] )
                            else: # 1d configuration
                                conc.append( composition[key] * np.ones(num_layers,dtype=np.float64) )

                        elif isinstance(composition[key], (list,np.ndarray)): # Already an array of initial conditions
                            conc.append( np.array(composition[key]),dtype=np.float64 )

                        elif isinstance(composition[key], dict): # Create array of initial conditions based off ratio to base element
                            factor_element = list(composition[key].keys())[0]   # This is the element that the initial condition is being based off of

                            if factor_element not in self.composition: # If the factor element is not yet in the composition list
                                self.composition.append(self.composition.pop(self.composition.index(key)))  # Bump current element to end of list
                                composition[key] = composition.pop(key)     # Bump current element to end of dictionary to pass through later

                            else: # Calculate the concentration of the current element
                                index = self.composition.index(factor_element)
                                if num_layers == 1: # 0d configuration
                                    conc.append( composition[key][factor_element] * conc[index] )
                                else: # 1d configuration
                                    conc.append( composition[key][factor_element] * conc[index] * np.ones(num_layers,dtype=np.float64) )

                else:
                    sys.exit("Phytoplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.conc = np.zeros((len(self.composition),num_layers-1,iters),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = scale * conc[const][:-1] # Apply scaling factor here to prevent from applying multiple times in the above step
        else:   # Model as single box
            self.conc = np.zeros((len(self.composition),iters),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = scale * conc[const]      # Apply scaling factor here to prevent from applying multiple times in the above step
        self.d_dt = np.zeros_like(self.conc[...,0],dtype=np.float64)
        self.conc_ratio = np.ones_like(self.conc[...,0],dtype=np.float64)

        # Add cell quotas ---------------------------------------------------------------
        # Create list of cell quota ids
        self.cell_quota_ids = List.empty_list(unicode_type)
        for element in self.composition:    # Base element not in cell quotas (cell quota of base element would be 1.)
            if element != base_element: self.cell_quota_ids.append(element)
        
        # Create list of maximum cell quotas
        if "max" in tracer["parameters"]["cell_quota"]: 
            self.cell_quota_max = List.empty_list(float64)
            for element in self.cell_quota_ids: # Add quotas in same order as ids
                self.cell_quota_max.append(tracer["parameters"]["cell_quota"]["max"][element])

        # Create list of minimum cell quotas
        if "min" in tracer["parameters"]["cell_quota"]: 
            self.cell_quota_min = List.empty_list(float64)
            for element in self.cell_quota_ids: # Add quotas in same order as ids
                self.cell_quota_min.append(tracer["parameters"]["cell_quota"]["min"][element])

        # Create list of optimal cell quotas
        if "opt" in tracer["parameters"]["cell_quota"]: 
            self.cell_quota_opt = List.empty_list(float64)
            for element in self.cell_quota_ids: # Add quotas in same order as ids
                self.cell_quota_opt.append(tracer["parameters"]["cell_quota"]["opt"][element])

        # Add production arrays ---------------------------------------------------------------
        self.upt = Dict.empty(
                key_type=types.unicode_type, 
                value_type=types.float64[:]
            )   # Uptake

        self.exu = np.zeros_like(self.conc[0,:,0],dtype=np.float64) # Exudation (Initialzied to 1 for use in respiration)
        self.gpp = np.zeros_like(self.conc[0,:,0],dtype=np.float64) # Gross Primary Production (Initialized to 1 for use in exudation and respiration)
        self.lys = np.zeros_like(self.conc[0,:,0],dtype=np.float64) # Lysis (carbon)
        self.npp = np.zeros_like(self.conc[0,...],dtype=np.float64) # Net Primary Production (include full time span for model output)
        self.psn = np.zeros_like(self.conc[0,:,0],dtype=np.float64) # Photosynthesis
        self.rsp = np.zeros_like(self.conc[0,:,0],dtype=np.float64) # Respiration

        # Add relevant reactions ---------------------------------------------------------------
        self.reactions = []
        for reac in reactions:
            # Add reaction to dictionary
            if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
            else:   consumed = {"empty": "empty"}
            if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
            else:   produced = {"empty": "empty"}
            if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
                self.reactions.append(reac)

        # Reorder uptake reactions in case of coupled uptake
        for i in range(len(self.reactions)):
            if self.reactions[i]["type"] == "uptake":
                # Get element of nutrient being consumed
                nutrient_element = list(self.reactions[i]["consumed"].values())[0][0]

                # Bump reaction to end of list if the nutrient is coupled to a different uptake rate
                if nutrient_element in self.coupled_uptake: self.reactions.append(self.reactions.pop(i))

        # Reorder reactions
        # uptake / gpp --> exu --> mortality --> respiration --> synthesis
        self.reactions = [item for item in self.reactions if item["type"] == "respiration"] + [item for item in self.reactions if item["type"] != "respiration"]
        self.reactions = [item for item in self.reactions if item["type"] == "mortality"] + [item for item in self.reactions if item["type"] != "mortality"]
        self.reactions = [item for item in self.reactions if item["type"] == "exudation"] + [item for item in self.reactions if item["type"] != "exudation"]
        self.reactions = [item for item in self.reactions if item["type"] == "gross_primary_production"] + [item for item in self.reactions if item["type"] != "gross_primary_production"]
        self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]
        self.reactions = [item for item in self.reactions if item["type"] == "photosynthesis"] + [item for item in self.reactions if item["type"] != "photosynthesis"]

        # Switch to determine if it is necessary to calculate growth parameters
        self.growth_switch = False
        for reac in self.reactions:
            if reac["type"] == "gross_primary_production" or reac["type"] == "uptake":
                self.growth_switch = True
                break

        # Boolean to determine whether activity and basal respiration will be calculated for chlorophyll synthesis
        self.calc_respiration = False
        for reac in self.reactions:
            if reac["type"] == "respiration":
                self.calc_respiration = True
                break


    def phyto(self, iter, base_element, temperature, z, dz, k_PAR, surface_PAR, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers, sinking):
        
        # Zero out variables
        for nut in self.upt:                                # Uptake
            self.upt[nut] = np.zeros_like(conc[0],dtype=np.float64)
        self.exu = np.zeros_like(conc[0],dtype=np.float64)  # Exudation (Initialzied to 1 for use in respiration)
        self.gpp = np.zeros_like(conc[0],dtype=np.float64)  # Gross Primary Production (Initialized to 1 for use in exudation and respiration)
        self.lys = np.zeros_like(conc[0],dtype=np.float64)  # Lysis (carbon)
        self.psn = np.zeros_like(conc[0],dtype=np.float64)  # Photosynthesis
        self.rsp = np.zeros_like(conc[0],dtype=np.float64)  # Respiration
        max_photo_rate = 0.

        # Initialize respiration arrays in case respiration is not included as a rate
        activity_respiration = np.zeros_like(conc[0],dtype=np.float64)
        basal_respiration = np.zeros_like(conc[0],dtype=np.float64)

        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(conc, conc_ratio, tracer_map, tracers)

        # Calculate temp regulation factor (if necessary)
        if self.temp_limited:   self.temp_regulation_factor = temperature_dependence(temperature, self.temp_reg_ids, self.temp_reg_params)
        
        # Calculate rates required for net primary production
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "exudation":                     self.exu = self.exudation(base_element, c, p, ec, ep, ic, ip, self.exudation_ids, self.exudation_params, self.nutrient_colimitation_factor, self.psn, self.upt, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "gross_primary_production":      self.gpp = self.gross_primary_production(self.abbrev, base_element, c, p, self.growth_ids, self.growth_params, self.psn, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "lysis":    
                composition_phyto = self.composition
                composition_om = tracers[p[0]].composition            
                self.lys = self.lysis(base_element, c, p, ec, ep, ic, ip, self.lysis_ids, self.lysis_params, self.lysis_apportioning_factor, self.cell_quota_ids, self.cell_quota_min, self.cell_quota_opt, self.nutrient_limitation["include"], self.nutrient_colimitation_factor, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, tracer_type, composition_phyto, composition_om)
            # if reac["type"] == "photosynthesis":                self.psn, fI, irr = self.photosynthesis(self.abbrev, self.growth_ids, self.growth_params, z, z, k_PAR, temperature, surface_PAR, self.temp_regulation_factor, self.nutrient_colimitation_factor, self.nutrient_limitation_factor, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "photosynthesis":                self.psn, irr, max_photo_rate = self.photosynthesis(self.abbrev, self.growth_ids, self.growth_params, z, z, k_PAR, temperature, surface_PAR, self.temp_regulation_factor, self.nutrient_colimitation_factor, self.nutrient_limitation_factor, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "respiration":                   self.rsp, activity_respiration, basal_respiration = self.respiration(self.abbrev, base_element, c, p, ec, ep, ic, ip, self.respiration_ids, self.respiration_params, self.temp_regulation_factor, self.exu, self.psn, conc, d_dt, tracer_map, self.composition)

        # Calculate net primary production
        base_index = self.composition.index(base_element)
        self.npp[:,iter] = self.net_primary_production(conc[tracer_map[self.abbrev][base_index]], self.exu, self.lys, self.psn, self.rsp)

        # Calculate remaining rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "chlorophyll_synthesis":         
                if self.calc_respiration:       self.chlorophyll_synthesis(self.abbrev, base_element, self.growth_ids, self.growth_params, activity_respiration, basal_respiration, irr, self.exu, self.lys, self.psn, conc, d_dt, tracer_map, self.composition)
                else:                           self.chlorophyll_synthesis(self.abbrev, base_element, self.growth_ids, self.growth_params, 0., 0., irr, self.exu, self.lys, self.psn, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "sedimentation": self.sedimentation(self.abbrev, self.background_sinking_rate, self.max_sinking_rate, self.sinking_threshold, self.nutrient_colimitation_factor, self.nutrient_limitation_factor, tracer_map, sinking, self.composition)
            if reac["type"] == "uptake": 
                # if c[0] == "no3" or c[0] == "nh4":  nh4_inhibited = self.nutrient_limitation["no3"]["nh4_inhibited"] 
                if c[0] == "no3" or c[0] == "nh4":  nh4_inhibited = self.nutrient_limitation[c[0]]["nh4_inhibited"] 
                else:   nh4_inhibited = False  
                # self.upt[c[0]] = self.uptake(self.abbrev, base_element, c, p, ec, ep, ic, self.uptake_ids, self.uptake_params, self.upt, self.coupled_uptake, self.cell_quota_ids, self.cell_quota_max, self.cell_quota_opt, self.temp_regulation_factor, self.nutrient_limitation_factor, nh4_inhibited, self.npp[...,iter], basal_respiration, self.psn, conc, conc_ratio, d_dt, tracer_map, self.composition)
                self.upt[c[0]] = self.uptake(self.abbrev, base_element, c, p, ec, ep, ic, self.uptake_ids, self.uptake_params, self.upt, self.coupled_uptake, self.cell_quota_ids, self.cell_quota_max, self.cell_quota_opt, self.temp_regulation_factor, self.nutrient_limitation_factor, nh4_inhibited, self.npp[...,iter], basal_respiration, self.psn, max_photo_rate, conc, conc_ratio, d_dt, tracer_map, self.composition)


    def add_nutrient(self, nutrients):
        """
        Add nutrients to phytoplankton and append dictionary of uptake rates
        """
        zeros = List.empty_list(float64[:])
        zeros.append(np.zeros(self.conc.shape[1],dtype=np.float64))
        for nut in nutrients:
            self.upt[nut] = np.zeros_like(self.conc[0,:,0],dtype=np.float64)
            # if nut in self.nutrient_limitation["include"]:  self.nutrient_limitation_factor[nut] = np.zeros_like(self.conc[0,:,0],dtype=np.float64)
            self.nutrient_limitation_factor[nut] = zeros


    def calculate_nutrient_limitation(self, conc, conc_ratio, tracer_map, tracers):
        """
        Definition:: Calculates nutrient limitation factor as either a minimum, product, or sum of all nutrients which limit phytoplankton growth
        """

        # fN = []
        fN = List.empty_list(float64[:])
        for key in self.nutrient_limitation:
            if key != "colimitation" and key != "include":
                # Get nutrient chemical constituent
                element = tracers[key].composition[0]

                # Get index of nutrient chemcical constituent in phytoplankton composition dictionary
                element_index = self.composition.index(element)

                # Get index of nutrient chemcical constituent in cell quota dictionary
                # quota_index = self.cell_quota.index(element)
                quota_index = self.cell_quota_ids.index(element)
                
                if self.nutrient_limitation[key]["type"] == "internal":
                    # Calculate nutrient limitation factor
                    func = ( conc_ratio[tracer_map[self.abbrev][element_index]] - self.cell_quota_min[quota_index]) / ( self.cell_quota_opt[quota_index] - self.cell_quota_min[quota_index] )
                    
                    # Ensures nonzero value
                    # func = np.maximum(1.E-20*np.ones_like(func), func)
                    func = np.minimum(np.ones_like(func),np.maximum(1.E-20*np.ones_like(func), func))   # maximum value of 1.

                    # Update dictionary
                    lim = List.empty_list(float64[:])
                    lim.append(np.minimum(np.ones_like(func),func))
                    self.nutrient_limitation_factor[key] = lim
                    # self.nutrient_limitation_factor[key] = [np.minimum(np.ones_like(func), func)]

                    # Append fN for colimitation calculation
                    # if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                    if key in self.nutrient_limitation["include"]:
                        fN.append(func)
                
                elif self.nutrient_limitation[key]["type"] == "external":
                    # Determine if Hill exponent exists for monod function
                    if "exponent" in self.nutrient_limitation[key]:
                        exponent = self.nutrient_limitation[key]["exponent"]
                    else:   # Default to 1. (no scaling)
                        exponent = 1.

                    # Calculate nutrient limitation factor
                    func = monod(conc[tracer_map[key][0]], self.nutrient_limitation[key]["half_sat"], exponent)
                    
                    # Ensures nonzero value
                    # func = np.maximum(1.E-20*np.ones_like(func), func)
                    func = np.minimum(np.ones_like(func),np.maximum(1.E-20*np.ones_like(func), func))   # maximum value of 1.

                    # Update dictionary
                    lim = List.empty_list(float64[:])
                    lim.append(np.minimum(np.ones_like(func),func))
                    self.nutrient_limitation_factor[key] = lim
                    # self.nutrient_limitation_factor[key] = np.minimum(np.ones_like(func), func)

                    # Append fN for colimitation calculation
                    # if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                    if "colimitation" in self.nutrient_limitation and key in self.nutrient_limitation["include"]:
                        fN.append(func)

                else:
                    sys.exit("Nutrient limitation type not recognized. Check documentation and edit input file.")

        if "colimitation" in self.nutrient_limitation.keys():   # Multiple nutrients available
            fN = np.array(fN)
            if self.nutrient_limitation["colimitation"] == "minimum":
                self.nutrient_colimitation_factor = np.min(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "product":
                self.nutrient_colimitation_factor = np.prod(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "sum":
                self.nutrient_colimitation_factor == np.sum(fN,  axis=0)
            else:
                sys.exit("Nutrient colimitation not recognized. Check documentation and edit input file.")
        else:   # Only one nutrient available
            self.nutrient_colimitation_factor = fN


    def aggregation(self, iter, base_element, parameters, c, ec, ic, tracers):
        """
        Definition:: Calculates the aggregation loss of phytoplankton during high biomass, bloom conditions
        Parameterized as a quadratic (density dependent) loss term
        """
        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Get concentration of base element
        index = self.composition.index(base_element)
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Calculate growth ratio
        growth_ratio = np.minimum(np.ones_like(self.psn[iter]), self.psn[iter]/ ( self.temp_regulation_factor * parameters["aggregation_frac"] * parameters["max_photo_rate"]))

        # Calculate aggregation limit
        agg_limit = (1. - growth_ratio)**2

        # Calculate aggregation loss
        agg_loss = agg_limit * parameters["aggegation_loss"] * (phyto**2)

        # Convert aggregation rate (if necessary)
        if "convert_aggregation" in parameters:
            if parameters["convert_aggregation"] == "cell_quota":
                quota_index = self.nutrient_limitation["nutrients"].index(c)
                agg_loss *= self.nutrient_limitation_factor[quota_index]
            else:
                if isinstance(parameters["convert_aggregation"],(int,float)) and not isinstance(parameters["convert_aggregation"],bool):
                    agg_loss *= parameters["convert_aggregation"]   
                elif isinstance(parameters["convert_aggregation"],str):   
                    agg_loss *= float(Fraction(parameters["convert_aggregation"]))

        # Calculate concentration ratio
        concentration_ratio(iter, index, self)

        # Update d_dt
        self.d_dt -= ec * tracers[c].conc_ratio * agg_loss

        # Apply partition to organic matter group (if necessary)
        if "partition" in parameters:   tracers[p].d_dt += ep * agg_loss * parameters["partition"]
        else:                           tracers[p].d_dt += ep * agg_loss


    @staticmethod
    @njit
    def chlorophyll_synthesis(abbrev, base_element, growth_ids, growth_params, activity_respiration, basal_respiration, irr, exudation, lysis, photosynthesis, conc, d_dt, tracer_map, composition):
        
        # Needs:    Chlorophyll quota (theta_chl in python bfm), Initial slope of PI curve (alpha_chl), Optimal value Epar/EK (p_EpEk_or),
        #           Maximal productivity (rp0), Chl:C relaxation rate (p_tochl_relt)

        # Extract parameter indices
        chl_quota = growth_ids.index("chl_quota")
        chl_relax_rate = growth_ids.index("chl_relax_rate")
        initial_PI_slope = growth_ids.index("initial_PI_slope")
        max_photo_rate = growth_ids.index("max_photo_rate")
        optimal_Epar_Ek = growth_ids.index("optimal_Epar_Ek")

        # Get concentration of base element
        base_index = composition.index(base_element)
        phyto = conc[tracer_map[abbrev][base_index]]

        # Get chlorophyll concentration
        chl_index = composition.index("chl")
        phyto_chl = conc[tracer_map[abbrev][chl_index]]

        # Calculate irradiance
        rho_chl = growth_params[chl_quota] * np.minimum(np.ones_like(photosynthesis), (photosynthesis - exudation - activity_respiration) * phyto / ( growth_params[initial_PI_slope] * ( phyto_chl + 1.E-20) * irr ))
        chl_opt = growth_params[optimal_Epar_Ek] * growth_params[max_photo_rate] * phyto / ( growth_params[initial_PI_slope] * irr + 1.E-20 )

        chlorophyll_synthesis = rho_chl * ( photosynthesis - exudation - activity_respiration ) * phyto \
                                    - ( lysis + basal_respiration ) * phyto_chl - np.maximum(np.zeros_like(phyto_chl), phyto_chl - chl_opt) * growth_params[chl_relax_rate]

        # Update d_dt
        d_dt[tracer_map[abbrev][chl_index]] += chlorophyll_synthesis


    @staticmethod
    @njit
    def exudation(base_element, c, p, ec, ep, ic, ip, exudation_ids, exudation_params, nutrient_colimitation_factor, photosynthesis, uptake, conc, d_dt, tracer_map, composition):

        # Extract parameter indices
        exc_frac = exudation_ids.index("excreted_fraction")
        method = exudation_ids.index("method")

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons][0]
        ind_p = ip[prod][0]

        # Concentration of consumed constituent
        tc = conc[tracer_map[cons][ind_c]]

        if exudation_params[method] == 1.: # "constant":
            exudation = exudation_params[exc_frac] * tc

            # Update d_dt
            d_dt[tracer_map[cons][ind_c]] -= exudation
            d_dt[tracer_map[prod][ind_p]] += exudation

            if composition[ind_c] == base_element:
                return exudation
                

        elif exudation_params[method] == 2.: # "photosynthesis":
            # Calculate activity and nutrient stress components
            activity = photosynthesis * exudation_params[exc_frac]
            # nutrient_stress = self.psn[iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_limitation_factor )
            nutrient_stress = photosynthesis * ( 1. - exudation_params[exc_frac] ) * ( 1. - nutrient_colimitation_factor )
            
            # Calculate exudation rate
            exu = activity + nutrient_stress
            exudation = exu * tc

            # Update d_dt
            d_dt[tracer_map[cons][ind_c]] -= exudation
            d_dt[tracer_map[prod][ind_p]] += exudation

            if composition[ind_c] == base_element:
                return exu

        elif exudation_params[method] == 3.: # "uptake":
            # Dictionary of potential nutrient compositions
            comps = Dict.empty(key_type=types.unicode_type,value_type=types.unicode_type)
            comps["no3"] = "n"
            comps["nh4"] = "n"
            comps["po4"] = "p"
            comps["fe"] = "fe"
            comps["sio4"] = "si"

            # Calculate total uptake rate for element
            uptake_sum = np.zeros_like(tc)
            for nut in uptake:
                if comps[nut] == composition[ind_c]:
                    uptake_sum += uptake[nut]
            
            # Calculate exudation rate
            exudation = exudation_params[exc_frac] * np.maximum(np.zeros_like(uptake_sum), uptake_sum)

            # Update d_dt
            d_dt[tracer_map[cons][ind_c]] -= exudation
            d_dt[tracer_map[prod][ind_p]] += exudation

            if composition[ind_c] == base_element:
                return exudation
    
        return np.zeros_like(tc)
    

    @staticmethod
    @njit
    def gross_primary_production(abbrev, base_element, c, p, growth_ids, growth_params, photosynthesis, conc, d_dt, tracer_map, composition):
        
        # Extract parameter indices
        if "o2" in p:
            conv_o2 = False
            if "convert_o2" in growth_ids:
                conv_o2 = True
                convert_o2 = growth_ids.index("convert_o2")

        if "co2" in c:
            conv_co2 = False
            if "convert_co2" in growth_ids:
                conv_co2 = True
                convert_co2 = growth_ids.index("convert_co2")

        # Locate index of base element
        base_index = composition.index(base_element)

        # Get base element concentration
        phyto = conc[tracer_map[abbrev][base_index]]

        # Calculate gross primary production
        gross_primary_production = photosynthesis * phyto

        # Update d_dt
        d_dt[tracer_map[abbrev][base_index]] += gross_primary_production

        if "o2" in p:   
            if conv_o2: d_dt[tracer_map["o2"][0]] += gross_primary_production * growth_params[convert_o2]
            else:       d_dt[tracer_map["o2"][0]] += gross_primary_production

        if "co2"in c:   
            if conv_co2:    d_dt[tracer_map["co2"][0]] -= gross_primary_production * growth_params[convert_co2]
            else:           d_dt[tracer_map["co2"][0]] -= gross_primary_production

        return gross_primary_production


    @staticmethod
    @njit
    def lysis(base_element, c, p, ec, ep, ic, ip, lysis_ids, lysis_params, lysis_apportioning_factor, cell_quota_ids, cell_quota_min, cell_quota_opt, nutrient_limitation, nutrient_colimitation_factor, temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, tracer_type, composition_phyto, composition_om):
        
        # Needs:    Nutrient limitation, Half sat for stress lysis (h_pnp from python bfm), Activity respiration fraction (d_P0 from python bfm), 
        #           Extra lysis rate (p_seo from python bfm), Half sat for extra lysis (p_sheo from python bfm)

        # Extract parameter indices
        method = lysis_ids.index("method")
        if lysis_params[method] == 1.:      # cell quota
            max_stress_lysis_rate = lysis_ids.index("max_stress_lysis_rate")
            extra_lysis_rate = lysis_ids.index("extra_lysis_rate")
            half_sat_stress_lysis = lysis_ids.index("half_sat_stress_lysis")
            half_sat_extra_lysis = lysis_ids.index("half_sat_extra_lysis")
        elif lysis_params[method] == 2.:    # constant lysis rate
            lysis_rate = lysis_ids.index("lysis_rate")
            convert = False
            part = False
            if "convert_lysis" in lysis_ids:
                convert = True
                convert_lysis = lysis_ids.index("convert_lysis")
            if "partition" in lysis_ids:
                part = True
                partition = lysis_ids.index("partition")
        
        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons]
        ind_p = ip[prod]

        # Locate index of base element
        base_index = composition_phyto.index(base_element)
        phyto = conc[tracer_map[cons][base_index]]

        if lysis_params[method] == 1.:  # "cell_quota"
            # Calculate element ratios
            if "n" in composition_phyto:
                nitrogen_index = composition_phyto.index("n")
                nit_base = conc[tracer_map[cons][nitrogen_index]] / phyto
            else:
                nit_base = np.ones_like(phyto)

            if "p" in composition_phyto: 
                phosphorus_index = composition_phyto.index("p")
                phos_base = conc[tracer_map[cons][phosphorus_index]] / phyto
            else:
                phos_base = np.ones_like(phyto)

            # Extract nutrient quotas
            if "no3" in nutrient_limitation:
                quota_index = cell_quota_ids.index("n")
                min_nitrogen_quota = cell_quota_min[quota_index]
            else:
                min_nitrogen_quota = 0.

            if "po4" in nutrient_limitation:
                quota_index = cell_quota_ids.index("p")
                min_phosphorus_quota = cell_quota_min[quota_index]
            else: 
                min_phosphorus_quota = 0.
            
            # Calculate fraction of lysis released to dissolved pool
            min_quota = np.minimum(min_nitrogen_quota/(nit_base + 1.E-20), min_phosphorus_quota/(phos_base + 1.E-20))
            apportioning_factor = np.minimum(np.ones_like(min_quota), min_quota)

            if tracer_type[tracer_map[prod][0]] == "dissolved":     apportioning_factor = 1 - np.minimum(np.ones_like(min_quota), min_quota)
            elif tracer_type[tracer_map[prod][0]] == "particulate": apportioning_factor = np.minimum(np.ones_like(min_quota), min_quota)

            # Calculate nutrient stress lysis
            nutrient_stress_lysis = ( lysis_params[max_stress_lysis_rate] * lysis_params[half_sat_stress_lysis] ) / ( nutrient_colimitation_factor + lysis_params[half_sat_stress_lysis] ) \
                                        + ( lysis_params[extra_lysis_rate] * phyto ) / ( phyto + lysis_params[half_sat_extra_lysis] + 1.E-20 )
            
            # # Apportion lysis between organic matter pools based on type == particulate or dissolved
            # if tracer_type[tracer_map[prod][0]] == "dissolved":     lysis = ( 1 - apportioning_factor ) * nutrient_stress_lysis * phyto
            # elif tracer_type[tracer_map[prod][0]] == "particulate": lysis = apportioning_factor * nutrient_stress_lysis * phyto

            # Update lysis variable
            lys = nutrient_stress_lysis

            # Match phyto and organic matter concentration ratios
            # ratios = np.zeros(len(tracers[p].conc_ratio))
            ratios = np.zeros((len(composition_om),conc.shape[1]))
            for const in composition_phyto:
                if const in composition_om:
                    index_phyto = composition_phyto.index(const)
                    index_om = composition_om.index(const)
                    ratios[index_om] = conc_ratio[tracer_map[cons][index_phyto]]

            # # Update d_dt
            # for i in range(len(elem_c)):
            #     d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * lysis
            # for j in range(len(elem_p)):
            #     d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * lysis


            # Apply apportioning factor if necessary
            if prod in lysis_apportioning_factor:
                necessary_constituents = lysis_apportioning_factor[prod]
                for i in range(len(elem_c)):
                    if composition_phyto[i] in necessary_constituents: d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * apportioning_factor * nutrient_stress_lysis * phyto
                    else:   d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * nutrient_stress_lysis * phyto
                for j in range(len(elem_p)):
                    if composition_om[j] in necessary_constituents: d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * apportioning_factor * nutrient_stress_lysis * phyto
                    else:   d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * nutrient_stress_lysis * phyto
            
            else:
                for i in range(len(elem_c)):    d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * nutrient_stress_lysis * phyto
                for j in range(len(elem_p)):    d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * nutrient_stress_lysis * phyto


        elif lysis_params[method] == 2.:    # "constant"
            lys = lysis_params[lysis_rate] * temp_regulation_factor * ( phyto**2 )

            mult = np.zeros((len(composition_phyto),conc.shape[1]))
            ratios = np.zeros((len(composition_om),conc.shape[1]))
            for const in composition_phyto:
                if const in composition_om:
                    index_phyto = composition_phyto.index(const)
                    index_om = composition_om.index(const)
                    ratios[index_om] = conc_ratio[tracer_map[cons][index_phyto]]

            # Convert lysis rate (if necessary)
            if convert:
                if lysis_params[convert_lysis] == -1:   # concentration ratio
                    pass    # mult is already concentration ratio
                else:   # constant
                    mult *= lysis_params[convert_lysis]
                    ratios *= lysis_params[convert_lysis]
                    # for const in composition_phyto:
                    #     mult[const] *= lysis_params[convert_lysis][const]
                    # for const in composition_om:
                    #     ratios[const] *= lysis_params[convert_lysis][const]

            if part:
                ratios *= lysis_params[partition]
                # for const in composition_om:
                #     ratios[const] *= lysis_params[partition][const]

            # Update d_dt
            for i in range(len(elem_c)):
                d_dt[tracer_map[cons][i]] -= elem_c[i] * mult[i] * lys
            for j in range(len(elem_p)):
                d_dt[tracer_map[prod][i]] += elem_p[j] * ratios[j] * lys

        return lys


    @staticmethod
    @njit
    def net_primary_production(phyto, exudation, lysis, photosynthesis, respiration):
            
        # Calculate losses
        specific_losses = exudation + respiration + lysis

        # Calculate net primary production
        npp = np.maximum( np.zeros_like(phyto), ( photosynthesis - specific_losses ) * phyto )

        return npp


    @staticmethod
    @njit
    def photosynthesis(abbrev, growth_ids, growth_params, coordinates, dz, k_PAR, temperature, surface_PAR, temp_regulation_factor, nutrient_colimitation_factor, nutrient_limitation_factor, conc, d_dt, tracer_map, composition):

        # Extract parameter indices
        max_photo_rate = growth_ids.index("max_photo_rate")
        light_lim = growth_ids.index("light_limitation")
        si_lim = growth_ids.index("silicate_limitation")

        # Maximal productivity
        if growth_params[max_photo_rate] == -1.:     # "eppley"
            Vm = max_growth_rate(growth_ids, growth_params, temperature)
        else:
            Vm = growth_params[max_photo_rate] * np.ones(len(temperature))
        
        # Light limitation
        if growth_params[light_lim] in [-1., -2., -3., -4.]:   # ["geider", "monod", "platt", "smith"]:
            eps_PAR = growth_ids.index("eps_PAR")

            # Calculate irradiance at surface
            irrad = irradiance(growth_params[eps_PAR], surface_PAR, coordinates, k_PAR)

            # Calculate light limitation
            irr, fI = light_limitation(abbrev, growth_ids, growth_params, dz, irrad, k_PAR, Vm, temp_regulation_factor, nutrient_colimitation_factor, conc, tracer_map, composition)
        else:
            fI = growth_params[light_lim] * np.ones(len(temperature))
            irr = 1.E-20 * np.ones(len(temperature))

        # Photosynthesis rate
        photosynthesis = temp_regulation_factor * Vm * fI

        if growth_params[si_lim] == 1.:    # Photosynthesis is limited by silicate
            photosynthesis *= nutrient_limitation_factor["sio4"][0]
        
        # return photosynthesis, fI, irr
        return photosynthesis, irr, Vm
    

    @staticmethod
    @njit
    def respiration(abbrev, base_element, c, p, ec, ep, ic, ip, respiration_ids, respiration_params, temp_regulation_factor, exudation, photosynthesis, conc, d_dt, tracer_map, composition):
        """
        Definition:: Calculates phytoplankton respiration
        """

        # Needs:    Activity and Basal respiration rates (gammaP in python bfm), Activity and Nutrient stress excretion

        # Extract parameter indices
        activity_respiration_frac = respiration_ids.index("activity_respiration_frac")
        basal_respiration_rate = respiration_ids.index("basal_respiration_rate")
        conv_o2 = False
        if "convert_o2" in respiration_ids:
            conv_o2 = True
            convert_o2 = respiration_ids.index("convert_o2")
        conv_co2 = False
        if "convert_co2" in respiration_ids:
            conv_co2 = True
            convert_co2 = respiration_ids.index("convert_co2")

        # Extract dict
        if p:   # "produced" not an empty list
            if len(p) > 1:  # if "co2" also included as produced element
                for index in range(len(p)):
                    if p[index] == "co2":   pass
                    else:   prod = p[index]
            else:   prod = p[0]
            elem_p = ep[prod]
            ind_p = ip[prod][0]

        # Locate index of base element
        base_index = composition.index(base_element)

        # Get base_element concentration
        phyto = conc[tracer_map[abbrev][base_index]]

        # Respiration
        activity_respiration = respiration_params[activity_respiration_frac] * ( photosynthesis - exudation )
        basal_respiration = temp_regulation_factor * respiration_params[basal_respiration_rate]
        total_respiration = activity_respiration + basal_respiration

        respiration = total_respiration * phyto

        # Update d_dt
        d_dt[tracer_map[abbrev][base_index]] -= respiration
        if "o2" in c:   
            if conv_o2:     d_dt[tracer_map["o2"][0]] -= respiration * respiration_params[convert_o2]
            else:           d_dt[tracer_map["o2"][0]] -= respiration

        if p:   # "produced" not an empty list
            for element in range(len(p)):
                if p[element] == "co2":
                    if base_element == "c":     d_dt[tracer_map["co2"][0]] += respiration  
                    else:                       
                        if conv_co2:    d_dt[tracer_map["co2"][0]] += respiration * respiration_params[convert_co2]
                        else:           d_dt[tracer_map["co2"][0]] += respiration
                else:   
                    d_dt[tracer_map[p[element]][ip[p[element]][0]]] += respiration

        return total_respiration, activity_respiration, basal_respiration


    @staticmethod
    @njit
    def sedimentation(abbrev, background_sinking_rate, max_sinking_rate, sinking_threshold, nutrient_colimitation_factor, nutrient_limitation_factor, tracer_map, sinking, composition):
        
        # Initialize sedimentation limiting factor
        sedi_lim = nutrient_colimitation_factor
        
        # Add silicate to sedi_lim if silicate is a phytoplankton constituent
        if "sio4" in tracer_map and "si" in composition:
            # Minimum between current colimiation factor and silicate limitation factor
            sedi_lim = np.minimum(sedi_lim, nutrient_limitation_factor["sio4"][0])

        for const in tracer_map[abbrev]:
            sinking[const] = background_sinking_rate + max_sinking_rate * np.maximum(np.zeros_like(sedi_lim), sinking_threshold - sedi_lim)
        
        return


    @staticmethod
    @njit
    def uptake(abbrev, base_element, c, p, ec, ep, ic, uptake_ids, uptake_params, upt, coupled_uptake_dict, cell_quota_ids, cell_quota_max, cell_quota_opt, temp_regulation_factor, nutrient_limitation_factor, nh4_inhibited, net_primary_production, basal_respiration, photosynthesis, max_photo_rate, conc, conc_ratio, d_dt, tracer_map, composition):
        
        element_compositions = Dict.empty(key_type=types.unicode_type,value_type=types.unicode_type)
        element_compositions["no3"] = "n"
        element_compositions["nh4"] = "n"
        element_compositions["po4"] = "p"
        element_compositions["fe"] = "fe"
        element_compositions["sio4"] = "si"

        # Identify the chemical constituent of the nutrient(s)
        element = c[0]

        # Identify uptake parameters for chemical constituent
        ids = uptake_ids[element]
        params = uptake_params[element]

        # Identify strategy for uptake rate calculation
        strategy = ids.index("strategy")

        # Get concentration of constituent in phytoplankton if present
        if element_compositions[element] in composition:     element_index = composition.index(element_compositions[element])

        if params[strategy] == 2.:  # "coupled" uptake
            # dictionary of codes for nutrient uptakes
            coupled_uptake_codes = Dict.empty(key_type=types.float64,value_type=types.unicode_type)
            coupled_uptake_codes = {np.float64(1.): "no3", np.float64(2.): "nh4", np.float64(3.): "po4", np.float64(4.): "fe", np.float64(5.): "sio4"}
            
            coupled_uptake = coupled_uptake_dict[element]
            linked_nutrients = coupled_uptake["links"]
            convert = False
            if "convert_uptake" in coupled_uptake:  convert = True

            # Extract uptake rates of linked nutrients
            uptake_rates = List.empty_list(float64[:])
            for nut in linked_nutrients:
                uptake_rates.append(upt[coupled_uptake_codes[nut]])

            # Calculate total linked uptake rate if multiple linked nutrients are used
            linked_uptake = uptake_rates[0].copy()
            if linked_nutrients and len(linked_nutrients) > 1:   # use numpy "maximum" to ensure minimum uptake of 0.
                if coupled_uptake["method"] == 1.:      # "max"
                    for i in range(1,len(uptake_rates)):
                        linked_uptake = np.maximum(linked_uptake,uptake_rates[i])

                elif coupled_uptake["method"] == 2.:    # "min"
                    for i in range(1,len(uptake_rates)):
                        linked_uptake = np.minimum(linked_uptake, uptake_rates[i])
                    
                elif coupled_uptake["method"] == 3.:    # "sum"
                    for i in range(1,len(uptake_rates)):
                        linked_uptake += uptake_rates[i]

                elif coupled_uptake["method"] == 4.:    # "product"
                    for i in range(1,len(uptake_rates)):
                        linked_uptake *= uptake_rates[i]
                
            # minimum uptake 0.
            linked_uptake = np.maximum(linked_uptake, np.zeros_like(uptake_rates[0]))

            if convert:
                convert_uptake = coupled_uptake["convert_uptake"][0]
                uptake = np.zeros(len(nutrient_limitation_factor[element][0]),dtype=np.float64)
                for depth in range(len(nutrient_limitation_factor[element][0])):
                    uptake[depth] = nutrient_limitation_factor[c[0]][0][depth] * linked_uptake[depth] * coupled_uptake["convert_uptake"][0]

            # Update d_dt
            for idx in range(len(c)):
                d_dt[tracer_map[c[idx]][0]] -= ec[c[idx]][0] * np.maximum(uptake, np.zeros_like(uptake))
            if element_compositions[element] in composition:     
                d_dt[tracer_map[abbrev][element_index]] += np.maximum(uptake, np.zeros_like(uptake))

        elif params[strategy] == 1.:  #"independent":
            
            cons = c[0]
            basis = ids.index("basis")

            if cons == "no3":
                # # Get concentration of element in phytoplankton
                if "n" in composition:
                    index = composition.index("n")
                else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                    index = composition.index(base_element)

                phyto = conc[tracer_map[abbrev][index]]

                # Get nutrient concentrations
                no3 = conc[tracer_map["no3"][0]]
                if "nh4" in tracer_map:     nh4 = conc[tracer_map["nh4"][0]]
                else:   nh4 = np.zeros_like(no3)

                # Determine uptake strategy
                if params[basis] == 1.:         # constant uptake rate
                    constant = ids.index("constant")
                    uptake = params[constant] * nutrient_limitation_factor["no3"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
                    uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["no3"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 3.:       # nutrient based uptake rate
                    if len(p) > 1: # If uptake can be source of organic matter
                        phy = list(p).index(abbrev)
                        ephy = ep[abbrev]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                        if phy == 0:    i = 1
                        else:           i = 0
                        om = p[i]
                        eom = ep[om]
                        om_nutrient_index = list(eom).index(1.)
                        om_nutrient = conc[tracer_map[om][om_nutrient_index]]
                    else:
                        phy = p[0]
                        ephy = ep[phy]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                    form = ids.index("form")
                    half_sat_nh4_preference = ids.index("half_sat_nh4_preference")
                    luxury_storage = ids.index("luxury_storage")
                    
                    # Get concentration of element in phytoplankton
                    if params[form] == 1.:  # "affinity":    # use specific affinity
                        specific_affinity = ids.index("specific_affinity")
                        index = composition.index(base_element)
                    else:   # use nutrient constituent
                        index = composition.index("n")

                    phyto = conc[tracer_map[abbrev][index]]

                    # Get nutrient concentrations
                    no3 = conc[tracer_map["no3"][0]]
                    if "nh4" in tracer_map:     nh4 = conc[tracer_map["nh4"][0]]
                    else:   nh4 = np.zeros_like(no3)

                    # Calculate preference for Ammonium uptake
                    nh4_preference = params[half_sat_nh4_preference] / ( params[half_sat_nh4_preference] + nh4 + 1.E-20)
                    
                    # Calculate maximum nitrogen uptake
                    max_uptake_no3 = params[specific_affinity] * no3 * phyto * nh4_preference
                    max_uptake_nh4 = params[specific_affinity] * nh4 * phyto
                    max_uptake_DIN = max_uptake_no3 + max_uptake_nh4

                    # Extract nutrient quota
                    # quota_index = self.nutrient_limitation["nutrients"].index("no3")
                    # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
                    quota_index = cell_quota_ids.index("n")
                    # nutrient_quota = cell_quota_opt[quota_index]
                    nutrient_quota = cell_quota_max[quota_index]

                    # Intracellular missing amount of N
                    missing_nit = max_photo_rate * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )
                    
                    # N uptake based on net assimilation of C
                    assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

                    # Actual uptake of nitrogen
                    actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

                    upt_switch = switch(actual_uptake)

                    no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
                    nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

                    # phyto_uptake = -actual_uptake * (1. - upt_switch)
                    uptake_to_phyto = no3_uptake + nh4_uptake
                    uptake_to_om = -actual_uptake * (1. - upt_switch)

                    uptake = no3_uptake
                    # if "nh4" in tracer_map:     uptake -= uptake_to_om/2    # /2 to split between no3 and nh4
                    # pass
                
                # Multiply by Monod function of nutrient limitation (if necessary)
                # if "nh4" in tracer_map and nh4_inhibited:    # no3_lim / (no3_lim + nh4_lim)
                #     uptake *= monod(nutrient_limitation_factor["no3"][0], nutrient_limitation_factor["nh4"][0], 1.)

                if "nh4" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
                    uptake_to_om /= 2

                # Update d_dt
                d_dt[tracer_map[cons][0]] -= uptake
                # if "n" in composition:     
                #     if "nh4" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
                #         d_dt[tracer_map[abbrev][index]] += uptake - uptake_to_om/2
                #     else:
                #         d_dt[tracer_map[abbrev][index]] += uptake - uptake_to_om
                    
                if len(p) > 1:  # both phyto and om produced
                    for i in range(len(ephy)):
                        d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)
                    for j in range(len(eom)):
                        d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                
                else:   # just phyto produced
                    for i in range(len(ephy)):
                        d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)

            elif cons == "nh4":
                # Get concentration of element in phytoplankton
                if "n" in composition:
                    index = composition.index("n")
                else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                    index = composition.index(base_element)

                phyto = conc[tracer_map[abbrev][index]]

                nh4 = conc[tracer_map["nh4"][0]]
                if "no3" in tracer_map:     no3 = conc[tracer_map["no3"][0]]
                else:   no3 = np.zeros_like(nh4)


                # Determine uptake strategy
                if params[basis] == 1.:         # constant uptake rate
                    constant = ids.index("constant")
                    uptake = params[constant] * nutrient_limitation_factor["nh4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
                    uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["nh4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 3.:       # nutrient based uptake rate
                    if len(p) > 1: # If uptake can be source of organic matter
                        phy = list(p).index(abbrev)
                        ephy = ep[abbrev]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                        if phy == 0:    i = 1
                        else:           i = 0
                        om = p[i]
                        eom = ep[om]
                        om_nutrient_index = list(eom).index(1.)
                        om_nutrient = conc[tracer_map[om][om_nutrient_index]]
                    else:
                        phy = p[0]
                        ephy = ep[phy]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                    form = ids.index("form")
                    half_sat_nh4_preference = ids.index("half_sat_nh4_preference")
                    luxury_storage = ids.index("luxury_storage")
                    
                    # Get concentration of element in phytoplankton
                    if params[form] == 1.:  # "affinity":    # use specific affinity
                        specific_affinity = ids.index("specific_affinity")
                        index = composition.index(base_element)
                    else:   # use nutrient constituent
                        index = composition.index("n")
                    
                    phyto = conc[tracer_map[abbrev][index]]

                    nh4 = conc[tracer_map["nh4"][0]]
                    if "no3" in tracer_map:     no3 = conc[tracer_map["no3"][0]]
                    else:   no3 = np.zeros_like(nh4)

                    # Calculate preference for Ammonium uptake
                    nh4_preference = params[half_sat_nh4_preference] / ( params[half_sat_nh4_preference] + nh4 + 1.E-20)

                    # Calculate maximum nitrogen uptake
                    max_uptake_no3 = params[specific_affinity] * no3 * phyto * nh4_preference
                    max_uptake_nh4 = params[specific_affinity] * nh4 * phyto
                    max_uptake_DIN = max_uptake_no3 + max_uptake_nh4
                    # Extract nutrient quota
                    # quota_index = self.nutrient_limitation["nutrients"].index("no3")
                    # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
                    quota_index = cell_quota_ids.index("n")
                    # nutrient_quota = cell_quota_opt[quota_index]
                    nutrient_quota = cell_quota_max[quota_index]

                    # Intracellular missing amount of N
                    missing_nit = max_photo_rate * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )
                    
                    # N uptake based on net assimilation of C
                    assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

                    # Actual uptake of nitrogen
                    actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

                    upt_switch = switch(actual_uptake)

                    no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
                    nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

                    # phyto_uptake = -actual_uptake * (1. - upt_switch)
                    uptake_to_phyto = no3_uptake + nh4_uptake
                    uptake_to_om = -actual_uptake * (1. - upt_switch)

                    uptake = nh4_uptake
                    # if "nh4" in tracer_map:     uptake -= uptake_to_om/2    # /2 to split between no3 and nh4
                    # pass

                # Multiply by Monod function of nutrient limitation (if necessary)
                # if nh4_inhibited:    # nh4_lim / (nh4_lim + no3_lim)
                #     uptake *= monod(nutrient_limitation_factor["nh4"][0], nutrient_limitation_factor["no3"][0], 1.)

                if "no3" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
                    uptake_to_om /= 2

                # Update d_dt
                d_dt[tracer_map[cons][0]] -= np.maximum(uptake, np.zeros_like(uptake))
                # if "n" in composition:     
                #     if "no3" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
                #         d_dt[tracer_map[abbrev][index]] += np.maximum(uptake, np.zeros_like(uptake)) - uptake_to_om/2
                #     else:
                #         d_dt[tracer_map[abbrev][index]] += np.maximum(uptake, np.zeros_like(uptake)) - uptake_to_om
                    
                if len(p) > 1:  # both phyto and om produced
                    for i in range(len(ephy)):
                        d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)
                    for j in range(len(eom)):
                        d_dt[tracer_map[om][j]] += eom[j] * uptake_to_om
                
                else:   # just phyto produced
                    for i in range(len(ephy)):
                        d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)

            elif cons == "po4":
                # Get concentration of element in phytoplankton
                # if "p" in self.composition:
                #     index = self.composition.index("p")
                # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                #     index = self.composition.index(base_element)
                form = ids.index("form")
                if params[form] ==  1.:     # use specific affinity
                    index = composition.index(base_element)
                else:   # use nutrient constituent
                    index = composition.index("p")
                
                phyto = conc[tracer_map[abbrev][index]]

                # Determine uptake strategy
                # Determine uptake strategy
                if params[basis] == 1.:         # constant uptake rate
                    constant = ids.index("constant")
                    uptake = params[constant] * nutrient_limitation_factor["po4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
                    uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["po4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)
                elif params[basis] == 3.:       # nutrient based uptake rate
                    if len(p) > 1: # If uptake can be source of organic matter
                        phy = list(p).index(abbrev)
                        ephy = ep[abbrev]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                        if phy == 0:    i = 1
                        else:           i = 0
                        om = p[i]
                        eom = ep[om]
                        om_nutrient_index = list(eom).index(1.)
                        om_nutrient = conc[tracer_map[om][om_nutrient_index]]
                    else:
                        phy = p[0]
                        ephy = ep[phy]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]
                    
                    # Get concentration of nutrient
                    # c = c[0]
                    elem_c = ec[cons][0]
                    ind_c = ic[cons][0]
                    nutrient = conc[tracer_map[cons][ind_c]]

                    # Calculate maximum nutrient uptake
                    specific_affinity = ids.index("specific_affinity")
                    max_uptake = params[specific_affinity] * nutrient * phyto

                    # Extract nutrient quota
                    quota_index = cell_quota_ids.index("p")
                    # nutrient_quota = cell_quota_opt[quota_index]
                    nutrient_quota = cell_quota_max[quota_index]

                    # Intracellular missing amount of nutrient
                    luxury_storage = ids.index("luxury_storage")
                    missing = max_photo_rate * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )

                    # Nutrient uptake based on net assimilation of C
                    assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

                    # Actual uptake of nutrient
                    actual_uptake = np.minimum(max_uptake, missing + assim_uptake)
                    # actual_uptake = np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))

                    upt_switch = switch(actual_uptake)

                    uptake_to_phyto = upt_switch * actual_uptake
                    uptake_to_om = -actual_uptake * (1. - upt_switch)

                    uptake = uptake_to_phyto

                # Update d_dt
                # tracers[c].d_dt -= np.array(ec) * uptake
                # tracers[self.abbrev].d_dt += ep * phyto_uptake
                d_dt[tracer_map[cons][0]] -= elem_c * uptake_to_phyto
                # tracers[self.abbrev].d_dt += ephy * uptake_to_phyto
                if len(p) > 1:  
                    # tracers[self.abbrev].d_dt -= ephy * uptake_to_om
                    # tracers[om].d_dt += eom * uptake_to_om
                    for i in range(len(ephy)):
                        # tracers[self.abbrev].d_dt[i] += ephy[i] * uptake_to_phyto
                        # tracers[self.abbrev].d_dt[i] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
                        d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
                        # d_dt[tracer_map[abbrev][i]] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
                        # tracers[self.abbrev].d_dt[i] += ephy[i] * np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))
                    for j in range(len(eom)):
                        d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                        
                else:
                    for i in range(len(ephy)):
                        d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
                        # d_dt[tracer_map[abbrev][i]] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
                        
                    

                    # if c == 'po4':  self.uptp[iter] = uptake_to_phyto
                
            elif cons == "fe":
                pass
                # convert = False
                # if "convert_uptake" in ids:     
                #     convert = True
                #     convert_uptake = ids.index("convert_uptake")

                # # Get concentration of element in phytoplankton
                # if "fe" in composition:
                #     index = composition.index("fe")
                # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                #     index = composition.index(base_element)

                # phyto = conc[tracer_map[abbrev][index]]

                # # Determine uptake strategy
                # if params[basis] == 1.:         # constant uptake rate
                #     constant = ids.index("constant")
                #     uptake = params[constant] * temp_regulation_factor * nutrient_limitation_factor["fe"] * phyto
                # elif params[basis] == 2.:       # uptake based on growth rate
                #     uptake = photosynthesis * nutrient_limitation_factor["fe"] * phyto
                # elif params[basis] == 3.:       # nutrient based uptake rate
                #     pass

                # # Mulitply by conversion factor (if needed)
                # if convert:     uptake *= params["convert_uptake"]

                # # Uptake is zero if Fe:Base ratio meets or exceeds maximum ratio
                # if "fe" in composition:    # Only need to calculate if iron is directly resolved
                #     iron_quota = cell_quota_ids.index("fe")
                #     if conc_ratio[tracer_map[abbrev][index]] >= cell_quota_max[iron_quota]:    uptake = np.zeros_like(uptake)

                # # Update d_dt
                # d_dt[tracer_map[c][0]] -= uptake
                # if "fe" in composition:    d_dt[tracer_map[abbrev][index]] += uptake

            elif cons == "sio4":    
                # Get concentration of element in phytoplankton
                # if "p" in self.composition:
                #     index = self.composition.index("p")
                # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                #     index = self.composition.index(base_element)
                form = ids.index("form")
                if params[form] ==  1.:     # use specific affinity
                    index = composition.index(base_element)
                else:   # use nutrient constituent
                    index = composition.index("si")
                
                phyto = conc[tracer_map[abbrev][index]]

                if len(p) > 1: # If uptake can be source of organic matter
                    phy = list(p).index(abbrev)
                    ephy = ep[abbrev]
                    if phy == 0:    i = 1
                    else:           i = 0
                    om = p[i]
                    eom = ep[om]
                else:
                    phy = p[0]
                    ephy = ep[phy]
                elem_c = ec[cons][0]
                    
                # Determine uptake strategy
                # Determine uptake strategy
                if params[basis] == 1.:         # constant uptake rate
                    constant = ids.index("constant")

                    if params[form] == 1.:  # specific affinity, use nutrient quota with concentration of base element
                        quota_index = cell_quota_ids.index("si")
                        nutrient_quota = cell_quota_opt[quota_index]
                        
                        uptake = params[constant] * nutrient_quota * phyto
                    else:   # use concentration of nutrient constituent
                        uptake = params[constant] * nutrient_limitation_factor["sio4"][0] * phyto
                    # uptake = params[constant] * nutrient_limitation_factor["sio4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)

                elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
                    if params[form] == 1.:  # specific affinity, use nutrient quota with concentration of base element
                        quota_index = cell_quota_ids.index("si")
                        nutrient_quota = cell_quota_opt[quota_index]
                        
                        uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_quota * phyto
                    else:   # use concentration of nutrient constituent
                        uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["sio4"][0] * phyto
                    # uptake = photosynthesis * nutrient_limitation_factor["sio4"][0] * phyto
                    uptake_to_om = np.zeros_like(uptake)

                elif params[basis] == 3.:       # nutrient based uptake rate
                    if len(p) > 1: # If uptake can be source of organic matter
                        phy = list(p).index(abbrev)
                        ephy = ep[abbrev]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

                        if phy == 0:    i = 1
                        else:           i = 0
                        om = p[i]
                        eom = ep[om]
                        om_nutrient_index = list(eom).index(1.)
                        om_nutrient = conc[tracer_map[om][om_nutrient_index]]
                    else:
                        phy = p[0]
                        ephy = ep[phy]
                        phyto_nutrient_index = list(ephy).index(1.)
                        phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]
                    
                    # Get concentration of nutrient
                    elem_c = ec[cons][0]
                    ind_c = ic[cons][0]
                    nutrient = conc[tracer_map[cons][ind_c]]

                    # Calculate maximum nutrient uptake
                    specific_affinity = ids.index("specific_affinity")
                    max_uptake = params[specific_affinity] * nutrient * phyto

                    # Extract nutrient quota
                    quota_index = cell_quota_ids.index("si")
                    # nutrient_quota = cell_quota_opt[quota_index]
                    nutrient_quota = cell_quota_max[quota_index]

                    # Intracellular missing amount of nutrient
                    luxury_storage = ids.index("luxury_storage")
                    missing = max_photo_rate * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )

                    # Nutrient uptake based on net assimilation of C
                    assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

                    # Actual uptake of nutrient
                    actual_uptake = np.minimum(max_uptake, missing + assim_uptake)

                    upt_switch = switch(actual_uptake)

                    uptake_to_phyto = upt_switch * actual_uptake
                    uptake_to_om = -actual_uptake * (1. - upt_switch)

                    uptake = uptake_to_phyto

                # Update d_dt
                d_dt[tracer_map[cons][0]] -= elem_c * uptake
                if len(p) > 1:  
                    for i in range(len(ephy)):  d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
                    for j in range(len(eom)):   d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                        
                else:
                    for i in range(len(ephy)):  d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)

        return uptake


    # @staticmethod
    # @njit
    # def uptake(abbrev, base_element, c, p, ec, ep, ic, uptake_ids, uptake_params, upt, coupled_uptake_dict, cell_quota_ids, cell_quota_max, cell_quota_opt, temp_regulation_factor, nutrient_limitation_factor, nh4_inhibited, net_primary_production, basal_respiration, photosynthesis, conc, conc_ratio, d_dt, tracer_map, composition):
        
    #     element_compositions = Dict.empty(key_type=types.unicode_type,value_type=types.unicode_type)
    #     element_compositions["no3"] = "n"
    #     element_compositions["nh4"] = "n"
    #     element_compositions["po4"] = "p"
    #     element_compositions["fe"] = "fe"
    #     element_compositions["sio4"] = "si"

    #     # Identify the chemical constituent of the nutrient(s)
    #     if c and len(c) > 1:    # if no3 and nh4 uptake happen in same rate, use no3 ids
    #         no3_index = c.index("no3")
    #         element = c[no3_index]
    #     else:   # use individual element
    #         element = c[0]

    #     # element = element_compositions[c[0]]
        
    #     # Identify uptake parameters for chemical constituent
    #     ids = uptake_ids[element]
    #     params = uptake_params[element]

    #     # Identify strategy for uptake rate calculation
    #     strategy = ids.index("strategy")

    #     # Get concentration of constituent in phytoplankton if present
    #     if element_compositions[element] in composition:     element_index = composition.index(element_compositions[element])

    #     if params[strategy] == 2.:  # "coupled" uptake
    #         # dictionary of codes for nutrient uptakes
    #         coupled_uptake_codes = Dict.empty(key_type=types.float64,value_type=types.unicode_type)
    #         coupled_uptake_codes = {np.float64(1.): "no3", np.float64(2.): "nh4", np.float64(3.): "po4", np.float64(4.): "fe", np.float64(5.): "sio4"}
            
    #         coupled_uptake = coupled_uptake_dict[element]
    #         linked_nutrients = coupled_uptake["links"]
    #         convert = False
    #         if "convert_uptake" in coupled_uptake:  convert = True

    #         # Extract uptake rates of linked nutrients
    #         uptake_rates = List.empty_list(float64[:])
    #         # uptake_rates = List()
    #         for nut in linked_nutrients:
    #             uptake_rates.append(upt[coupled_uptake_codes[nut]])

    #         # Calculate total linked uptake rate if multiple linked nutrients are used
    #         linked_uptake = uptake_rates[0].copy()
    #         if linked_nutrients and len(linked_nutrients) > 1:   # use numpy "maximum" to ensure minimum uptake of 0.
    #             if coupled_uptake["method"] == 1.:      # "max":       
    #                 # linked_uptake = np.maximum(np.maximum(uptake_rates), np.zeros_like(uptake_rates[0]))
    #                 # linked_uptake = uptake_rates[0].copy()
    #                 for i in range(1,len(uptake_rates)):
    #                     linked_uptake = np.maximum(linked_uptake,uptake_rates[i])

    #             elif coupled_uptake["method"] == 2.:    # "min":     
    #                 # linked_uptake = np.maximum(np.minimum(uptake_rates), np.zeros_like(uptake_rates[0]))
    #                 # linked_uptake = uptake_rates[0].copy()
    #                 for i in range(1,len(uptake_rates)):
    #                     linked_uptake = np.minimum(linked_uptake, uptake_rates[i])
                    
    #             elif coupled_uptake["method"] == 3.:    # "sum":     
    #                 # linked_uptake = np.maximum(np.sum(uptake_rates), np.zeros_like(uptake_rates[0]))
    #                 # linked_uptake = uptake_rates[0].copy()
    #                 for i in range(1,len(uptake_rates)):
    #                     linked_uptake += uptake_rates[i]

    #             elif coupled_uptake["method"] == 4.:    # "product": 
    #                 # linked_uptake = np.maximum(np.prod(uptake_rates), np.zeros_like(uptake_rates[0]))
    #                 # linked_uptake = uptake_rates[0].copy()
    #                 for i in range(1,len(uptake_rates)):
    #                     linked_uptake *= uptake_rates[i]
                
    #         # minimum uptake 0.
    #         linked_uptake = np.maximum(linked_uptake, np.zeros_like(uptake_rates[0]))

    #         if convert:
    #             convert_uptake = coupled_uptake["convert_uptake"][0]
    #             uptake = np.zeros(len(nutrient_limitation_factor[element][0]),dtype=np.float64)
    #             # if element == "n":  
    #             #     for depth in range(len(nutrient_limitation_factor[c[0]])):
    #             #         uptake[depth] = nutrient_limitation_factor["no3"][0][depth] * linked_uptake[depth] * convert_uptake # coupled_uptake["convert_uptake"][0]
    #             # else:   
    #             #     for depth in range(len(nutrient_limitation_factor[c[0]])):
    #             #         uptake[depth] = nutrient_limitation_factor[c[0]][0][depth] * linked_uptake[depth] * coupled_uptake["convert_uptake"]
                
    #             for depth in range(len(nutrient_limitation_factor[element][0])):
    #                 uptake[depth] = nutrient_limitation_factor[c[0]][0][depth] * linked_uptake[depth] * coupled_uptake["convert_uptake"][0]

    #         # Update d_dt
    #         for idx in range(len(c)):
    #             d_dt[tracer_map[c[idx]][0]] -= ec[c[idx]][0] * np.maximum(uptake, np.zeros_like(uptake))
    #         if element_compositions[element] in composition:     
    #             d_dt[tracer_map[abbrev][element_index]] += np.maximum(uptake, np.zeros_like(uptake))

    #     elif params[strategy] == 1.:  #"independent":
            
    #         if len(c) > 1:  # Used if no3 and nh4 are consumed together rather than individually
    #             if len(p) > 1:  # If uptake can be source of organic matter
    #                 phy = list(p).index(abbrev)
    #                 ephy = ep[abbrev]
    #                 phyto_nutrient_index = list(ephy).index(1.)
    #                 phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                 if phy == 0:    i = 1
    #                 else:           i = 0
    #                 om = p[i]
    #                 eom = ep[om]
    #                 om_nutrient_index = list(eom).index(1.)
    #                 om_nutrient = conc[tracer_map[abbrev][om_nutrient_index]]

    #             else:
    #                 phy = p[0]
    #                 ephy = ep[phy]
    #                 phyto_nutrient_index = list(ephy).index(1.)
    #                 phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #             form = ids.index("form")
    #             half_sat_nh4_preference = ids.index("half_sat_nh4_preference")
    #             luxury_storage = ids.index("luxury_storage")
    #             max_photo_rate = ids.index("max_photo_rate")

    #             # Get concentration of element in phytoplankton
    #             if params[form] == 1.:  # "affinity":    # use specific affinity
    #                 specific_affinity = ids.index("specific_affinity")
    #                 index = composition.index(base_element)
    #             else:   # use nutrient constituent
    #                 index = composition.index("n")
    #             # if "n" in self.composition:
    #             #     index = self.composition.index("n")
    #             # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #             #     index = self.composition.index(base_element)
    #             phyto = conc[tracer_map[abbrev][index]]

    #             # Get nutrient concentrations
    #             no3 = conc[tracer_map["no3"][0]]
    #             nh4 = conc[tracer_map["nh4"][0]]
                
    #             # Calculate preference for Ammonium uptake
    #             nh4_preference = params[half_sat_nh4_preference] / ( params[half_sat_nh4_preference] + nh4 + 1.E-20)

    #             # Calculate maximum nitrogen uptake
    #             max_uptake_no3 = params[specific_affinity] * no3 * phyto * nh4_preference
    #             max_uptake_nh4 = params[specific_affinity] * nh4 * phyto
    #             max_uptake_DIN = max_uptake_no3 + max_uptake_nh4

    #             # Extract nutrient quota
    #             # quota_index = self.nutrient_limitation["nutrients"].index("no3")
    #             # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
    #             quota_index = cell_quota_ids.index("n")
    #             # nutrient_quota = cell_quota_opt[quota_index]
    #             nutrient_quota = cell_quota_max[quota_index]

    #             # Intracellular missing amount of N
    #             missing_nit = params[max_photo_rate] * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )
                
    #             # N uptake based on net assimilation of C
    #             assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

    #             # Actual uptake of nitrogen
    #             actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

    #             upt_switch = switch(actual_uptake)

    #             no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
    #             nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

    #             # phyto_uptake = -actual_uptake * (1. - upt_switch)
    #             uptake_to_phyto = no3_uptake + nh4_uptake
    #             uptake_to_om = -actual_uptake * (1. - upt_switch)

    #             # Calculate n2 uptake for nitrogen fixers (if necessary)
    #             if "n2" in c:
    #                 n2_uptake = np.zeros_like(phyto_nutrient)
    #                 n2_uptake = ( 1. - nutrient_limitation_factor["no3"][0] - nutrient_limitation_factor["nh4"][0] ) * photosynthesis * phyto_nutrient
    #                 uptake_to_phyto += np.maximum(n2_uptake, np.zeros_like(n2_uptake))

    #             # Update d_dt
    #             d_dt[tracer_map["no3"][0]] -= no3_uptake
    #             d_dt[tracer_map["nh4"][0]] -= nh4_uptake
    #             d_dt[tracer_map[abbrev][phyto_nutrient_index]] += uptake_to_phyto
    #             if len(p) > 1:  
    #                 # this ephy line might be unnecessary so i'm commenting it out but leaving just in case
    #                 for i in range(len(ephy)):
    #                     d_dt[tracer_map[abbrev][i]] -= ephy[i] * uptake_to_om
    #                 for j in range(len(eom)):
    #                     d_dt[tracer_map[om][j]] += eom[j] * uptake_to_om

    #         else:
    #             cons = c[0]
    #             basis = ids.index("basis")

    #             if cons == "no3":
    #                 # # Get concentration of element in phytoplankton
    #                 # if "n" in composition:
    #                 #     index = composition.index("n")
    #                 # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #                 #     index = composition.index(base_element)

    #                 # phyto = conc[tracer_map[abbrev][index]]

    #                 # # Get nutrient concentrations
    #                 # no3 = conc[tracer_map["no3"][0]]
    #                 # if "nh4" in tracer_map:     nh4 = conc[tracer_map["nh4"][0]]
    #                 # else:   nh4 = np.zeros_like(no3)

    #                 # Determine uptake strategy
    #                 if params[basis] == 1.:         # constant uptake rate
    #                     constant = ids.index("constant")
    #                     uptake = params[constant] * nutrient_limitation_factor["no3"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
    #                     uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["no3"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 3.:       # nutrient based uptake rate
    #                     if len(p) > 1: # If uptake can be source of organic matter
    #                         phy = list(p).index(abbrev)
    #                         ephy = ep[abbrev]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                         if phy == 0:    i = 1
    #                         else:           i = 0
    #                         om = p[i]
    #                         eom = ep[om]
    #                         om_nutrient_index = list(eom).index(1.)
    #                         om_nutrient = conc[tracer_map[om][om_nutrient_index]]
    #                     else:
    #                         phy = p[0]
    #                         ephy = ep[phy]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                     form = ids.index("form")
    #                     half_sat_nh4_preference = ids.index("half_sat_nh4_preference")
    #                     luxury_storage = ids.index("luxury_storage")
    #                     max_photo_rate = ids.index("max_photo_rate")
                        
    #                     # Get concentration of element in phytoplankton
    #                     if params[form] == 1.:  # "affinity":    # use specific affinity
    #                         specific_affinity = ids.index("specific_affinity")
    #                         index = composition.index(base_element)
    #                     else:   # use nutrient constituent
    #                         index = composition.index("n")

    #                     phyto = conc[tracer_map[abbrev][index]]

    #                     # Get nutrient concentrations
    #                     no3 = conc[tracer_map["no3"][0]]
    #                     if "nh4" in tracer_map:     nh4 = conc[tracer_map["nh4"][0]]
    #                     else:   nh4 = np.zeros_like(no3)

    #                     # Calculate preference for Ammonium uptake
    #                     nh4_preference = params[half_sat_nh4_preference] / ( params[half_sat_nh4_preference] + nh4 + 1.E-20)
                        
    #                     # Calculate maximum nitrogen uptake
    #                     max_uptake_no3 = params[specific_affinity] * no3 * phyto * nh4_preference
    #                     max_uptake_nh4 = params[specific_affinity] * nh4 * phyto
    #                     max_uptake_DIN = max_uptake_no3 + max_uptake_nh4

    #                     # Extract nutrient quota
    #                     # quota_index = self.nutrient_limitation["nutrients"].index("no3")
    #                     # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
    #                     quota_index = cell_quota_ids.index("n")
    #                     # nutrient_quota = cell_quota_opt[quota_index]
    #                     nutrient_quota = cell_quota_max[quota_index]

    #                     # Intracellular missing amount of N
    #                     missing_nit = params[max_photo_rate] * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )
                        
    #                     # N uptake based on net assimilation of C
    #                     assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

    #                     # Actual uptake of nitrogen
    #                     actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

    #                     upt_switch = switch(actual_uptake)

    #                     no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
    #                     nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

    #                     # phyto_uptake = -actual_uptake * (1. - upt_switch)
    #                     uptake_to_phyto = no3_uptake + nh4_uptake
    #                     uptake_to_om = -actual_uptake * (1. - upt_switch)

    #                     uptake = no3_uptake
    #                     # if "nh4" in tracer_map:     uptake -= uptake_to_om/2    # /2 to split between no3 and nh4
    #                     # pass
                    
    #                 # Multiply by Monod function of nutrient limitation (if necessary)
    #                 # if "nh4" in tracer_map and nh4_inhibited:    # no3_lim / (no3_lim + nh4_lim)
    #                 #     uptake *= monod(nutrient_limitation_factor["no3"][0], nutrient_limitation_factor["nh4"][0], 1.)

    #                 if "nh4" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
    #                     uptake_to_om /= 2

    #                 # Update d_dt
    #                 d_dt[tracer_map[cons][0]] -= uptake
    #                 # if "n" in composition:     
    #                 #     if "nh4" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
    #                 #         d_dt[tracer_map[abbrev][index]] += uptake - uptake_to_om/2
    #                 #     else:
    #                 #         d_dt[tracer_map[abbrev][index]] += uptake - uptake_to_om
                        
    #                 if len(p) > 1:  # both phyto and om produced
    #                     for i in range(len(ephy)):
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)
    #                     for j in range(len(eom)):
    #                         d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                    
    #                 else:   # just phyto produced
    #                     for i in range(len(ephy)):
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)

    #             elif cons == "nh4":
    #                 # # Get concentration of element in phytoplankton
    #                 # if "n" in composition:
    #                 #     index = composition.index("n")
    #                 # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #                 #     index = composition.index(base_element)

    #                 # phyto = conc[tracer_map[abbrev][index]]

    #                 # nh4 = conc[tracer_map["nh4"][0]]
    #                 # if "no3" in tracer_map:     no3 = conc[tracer_map["no3"][0]]
    #                 # else:   no3 = np.zeros_like(nh4)


    #                 # Determine uptake strategy
    #                 if params[basis] == 1.:         # constant uptake rate
    #                     constant = ids.index("constant")
    #                     uptake = params[constant] * nutrient_limitation_factor["nh4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
    #                     uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["nh4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 3.:       # nutrient based uptake rate
    #                     if len(p) > 1: # If uptake can be source of organic matter
    #                         phy = list(p).index(abbrev)
    #                         ephy = ep[abbrev]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                         if phy == 0:    i = 1
    #                         else:           i = 0
    #                         om = p[i]
    #                         eom = ep[om]
    #                         om_nutrient_index = list(eom).index(1.)
    #                         om_nutrient = conc[tracer_map[om][om_nutrient_index]]
    #                     else:
    #                         phy = p[0]
    #                         ephy = ep[phy]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                     form = ids.index("form")
    #                     half_sat_nh4_preference = ids.index("half_sat_nh4_preference")
    #                     luxury_storage = ids.index("luxury_storage")
    #                     max_photo_rate = ids.index("max_photo_rate")
                        
    #                     # Get concentration of element in phytoplankton
    #                     if params[form] == 1.:  # "affinity":    # use specific affinity
    #                         specific_affinity = ids.index("specific_affinity")
    #                         index = composition.index(base_element)
    #                     else:   # use nutrient constituent
    #                         index = composition.index("n")
                        
    #                     phyto = conc[tracer_map[abbrev][index]]

    #                     nh4 = conc[tracer_map["nh4"][0]]
    #                     if "no3" in tracer_map:     no3 = conc[tracer_map["no3"][0]]
    #                     else:   no3 = np.zeros_like(nh4)

    #                     # Calculate preference for Ammonium uptake
    #                     nh4_preference = params[half_sat_nh4_preference] / ( params[half_sat_nh4_preference] + nh4 + 1.E-20)

    #                     # Calculate maximum nitrogen uptake
    #                     max_uptake_no3 = params[specific_affinity] * no3 * phyto * nh4_preference
    #                     max_uptake_nh4 = params[specific_affinity] * nh4 * phyto
    #                     max_uptake_DIN = max_uptake_no3 + max_uptake_nh4
    #                     # Extract nutrient quota
    #                     # quota_index = self.nutrient_limitation["nutrients"].index("no3")
    #                     # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
    #                     quota_index = cell_quota_ids.index("n")
    #                     # nutrient_quota = cell_quota_opt[quota_index]
    #                     nutrient_quota = cell_quota_max[quota_index]

    #                     # Intracellular missing amount of N
    #                     missing_nit = params[max_photo_rate] * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )
                        
    #                     # N uptake based on net assimilation of C
    #                     assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

    #                     # Actual uptake of nitrogen
    #                     actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

    #                     upt_switch = switch(actual_uptake)

    #                     no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
    #                     nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

    #                     # phyto_uptake = -actual_uptake * (1. - upt_switch)
    #                     uptake_to_phyto = no3_uptake + nh4_uptake
    #                     uptake_to_om = -actual_uptake * (1. - upt_switch)

    #                     uptake = nh4_uptake
    #                     # if "nh4" in tracer_map:     uptake -= uptake_to_om/2    # /2 to split between no3 and nh4
    #                     # pass

    #                 # Multiply by Monod function of nutrient limitation (if necessary)
    #                 # if nh4_inhibited:    # nh4_lim / (nh4_lim + no3_lim)
    #                 #     uptake *= monod(nutrient_limitation_factor["nh4"][0], nutrient_limitation_factor["no3"][0], 1.)

    #                 if "no3" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
    #                     uptake_to_om /= 2

    #                 # Update d_dt
    #                 d_dt[tracer_map[cons][0]] -= np.maximum(uptake, np.zeros_like(uptake))
    #                 # if "n" in composition:     
    #                 #     if "no3" in tracer_map: # /2 to split uptake_to_om between no3 and nh4 uptake rates
    #                 #         d_dt[tracer_map[abbrev][index]] += np.maximum(uptake, np.zeros_like(uptake)) - uptake_to_om/2
    #                 #     else:
    #                 #         d_dt[tracer_map[abbrev][index]] += np.maximum(uptake, np.zeros_like(uptake)) - uptake_to_om
                        
    #                 if len(p) > 1:  # both phyto and om produced
    #                     for i in range(len(ephy)):
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)
    #                     for j in range(len(eom)):
    #                         d_dt[tracer_map[om][j]] += eom[j] * uptake_to_om
                    
    #                 else:   # just phyto produced
    #                     for i in range(len(ephy)):
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] * (uptake - uptake_to_om)

    #             elif cons == "po4":
    #                 # Get concentration of element in phytoplankton
    #                 # if "p" in self.composition:
    #                 #     index = self.composition.index("p")
    #                 # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #                 #     index = self.composition.index(base_element)
    #                 form = ids.index("form")
    #                 if params[form] ==  1.:     # use specific affinity
    #                     index = composition.index(base_element)
    #                 else:   # use nutrient constituent
    #                     index = composition.index("p")
                    
    #                 phyto = conc[tracer_map[abbrev][index]]

    #                 # Determine uptake strategy
    #                 # Determine uptake strategy
    #                 if params[basis] == 1.:         # constant uptake rate
    #                     constant = ids.index("constant")
    #                     uptake = params[constant] * nutrient_limitation_factor["po4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
    #                     uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["po4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)
    #                 elif params[basis] == 3.:       # nutrient based uptake rate
    #                     if len(p) > 1: # If uptake can be source of organic matter
    #                         phy = list(p).index(abbrev)
    #                         ephy = ep[abbrev]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                         if phy == 0:    i = 1
    #                         else:           i = 0
    #                         om = p[i]
    #                         eom = ep[om]
    #                         om_nutrient_index = list(eom).index(1.)
    #                         om_nutrient = conc[tracer_map[om][om_nutrient_index]]
    #                     else:
    #                         phy = p[0]
    #                         ephy = ep[phy]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]
                        
    #                     # Get concentration of nutrient
    #                     # c = c[0]
    #                     elem_c = ec[cons][0]
    #                     ind_c = ic[cons][0]
    #                     nutrient = conc[tracer_map[cons][ind_c]]

    #                     # Calculate maximum nutrient uptake
    #                     specific_affinity = ids.index("specific_affinity")
    #                     max_uptake = params[specific_affinity] * nutrient * phyto

    #                     # Extract nutrient quota
    #                     quota_index = cell_quota_ids.index("p")
    #                     # nutrient_quota = cell_quota_opt[quota_index]
    #                     nutrient_quota = cell_quota_max[quota_index]

    #                     # Intracellular missing amount of nutrient
    #                     max_photo_rate = ids.index("max_photo_rate")
    #                     luxury_storage = ids.index("luxury_storage")
    #                     missing = params[max_photo_rate] * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )

    #                     # Nutrient uptake based on net assimilation of C
    #                     assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

    #                     # Actual uptake of nutrient
    #                     actual_uptake = np.minimum(max_uptake, missing + assim_uptake)
    #                     # actual_uptake = np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))

    #                     upt_switch = switch(actual_uptake)

    #                     uptake_to_phyto = upt_switch * actual_uptake
    #                     uptake_to_om = -actual_uptake * (1. - upt_switch)

    #                     uptake = uptake_to_phyto

    #                 # Update d_dt
    #                 # tracers[c].d_dt -= np.array(ec) * uptake
    #                 # tracers[self.abbrev].d_dt += ep * phyto_uptake
    #                 d_dt[tracer_map[cons][0]] -= elem_c * uptake_to_phyto
    #                 # tracers[self.abbrev].d_dt += ephy * uptake_to_phyto
    #                 if len(p) > 1:  
    #                     # tracers[self.abbrev].d_dt -= ephy * uptake_to_om
    #                     # tracers[om].d_dt += eom * uptake_to_om
    #                     for i in range(len(ephy)):
    #                         # tracers[self.abbrev].d_dt[i] += ephy[i] * uptake_to_phyto
    #                         # tracers[self.abbrev].d_dt[i] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
    #                         # d_dt[tracer_map[abbrev][i]] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
    #                         # tracers[self.abbrev].d_dt[i] += ephy[i] * np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))
    #                     for j in range(len(eom)):
    #                         d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                            
    #                 else:
    #                     for i in range(len(ephy)):
    #                         d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
    #                         # d_dt[tracer_map[abbrev][i]] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
                            
                        

    #                     # if c == 'po4':  self.uptp[iter] = uptake_to_phyto
                    
    #             elif cons == "fe":
    #                 pass
    #                 # convert = False
    #                 # if "convert_uptake" in ids:     
    #                 #     convert = True
    #                 #     convert_uptake = ids.index("convert_uptake")

    #                 # # Get concentration of element in phytoplankton
    #                 # if "fe" in composition:
    #                 #     index = composition.index("fe")
    #                 # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #                 #     index = composition.index(base_element)

    #                 # phyto = conc[tracer_map[abbrev][index]]

    #                 # # Determine uptake strategy
    #                 # if params[basis] == 1.:         # constant uptake rate
    #                 #     constant = ids.index("constant")
    #                 #     uptake = params[constant] * temp_regulation_factor * nutrient_limitation_factor["fe"] * phyto
    #                 # elif params[basis] == 2.:       # uptake based on growth rate
    #                 #     uptake = photosynthesis * nutrient_limitation_factor["fe"] * phyto
    #                 # elif params[basis] == 3.:       # nutrient based uptake rate
    #                 #     pass

    #                 # # Mulitply by conversion factor (if needed)
    #                 # if convert:     uptake *= params["convert_uptake"]

    #                 # # Uptake is zero if Fe:Base ratio meets or exceeds maximum ratio
    #                 # if "fe" in composition:    # Only need to calculate if iron is directly resolved
    #                 #     iron_quota = cell_quota_ids.index("fe")
    #                 #     if conc_ratio[tracer_map[abbrev][index]] >= cell_quota_max[iron_quota]:    uptake = np.zeros_like(uptake)

    #                 # # Update d_dt
    #                 # d_dt[tracer_map[c][0]] -= uptake
    #                 # if "fe" in composition:    d_dt[tracer_map[abbrev][index]] += uptake

    #             elif cons == "sio4":    
    #                 # Get concentration of element in phytoplankton
    #                 # if "p" in self.composition:
    #                 #     index = self.composition.index("p")
    #                 # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
    #                 #     index = self.composition.index(base_element)
    #                 form = ids.index("form")
    #                 if params[form] ==  1.:     # use specific affinity
    #                     index = composition.index(base_element)
    #                 else:   # use nutrient constituent
    #                     index = composition.index("si")
                    
    #                 phyto = conc[tracer_map[abbrev][index]]

    #                 if len(p) > 1: # If uptake can be source of organic matter
    #                     phy = list(p).index(abbrev)
    #                     ephy = ep[abbrev]
    #                     if phy == 0:    i = 1
    #                     else:           i = 0
    #                     om = p[i]
    #                     eom = ep[om]
    #                 else:
    #                     phy = p[0]
    #                     ephy = ep[phy]
    #                 elem_c = ec[cons][0]
                        
    #                 # Determine uptake strategy
    #                 # Determine uptake strategy
    #                 if params[basis] == 1.:         # constant uptake rate
    #                     constant = ids.index("constant")

    #                     if params[form] == 1.:  # specific affinity, use nutrient quota with concentration of base element
    #                         quota_index = cell_quota_ids.index("si")
    #                         nutrient_quota = cell_quota_opt[quota_index]
                            
    #                         uptake = params[constant] * nutrient_quota * phyto
    #                     else:   # use concentration of nutrient constituent
    #                         uptake = params[constant] * nutrient_limitation_factor["sio4"][0] * phyto
    #                     # uptake = params[constant] * nutrient_limitation_factor["sio4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)

    #                 elif params[basis] == 2.:       # uptake based on growth rate (excluding respiratory costs)
    #                     if params[form] == 1.:  # specific affinity, use nutrient quota with concentration of base element
    #                         quota_index = cell_quota_ids.index("si")
    #                         nutrient_quota = cell_quota_opt[quota_index]
                            
    #                         uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_quota * phyto
    #                     else:   # use concentration of nutrient constituent
    #                         uptake = np.maximum(np.zeros_like(photosynthesis,dtype=np.float64), photosynthesis - basal_respiration) * nutrient_limitation_factor["sio4"][0] * phyto
    #                     # uptake = photosynthesis * nutrient_limitation_factor["sio4"][0] * phyto
    #                     uptake_to_om = np.zeros_like(uptake)

    #                 elif params[basis] == 3.:       # nutrient based uptake rate
    #                     if len(p) > 1: # If uptake can be source of organic matter
    #                         phy = list(p).index(abbrev)
    #                         ephy = ep[abbrev]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]

    #                         if phy == 0:    i = 1
    #                         else:           i = 0
    #                         om = p[i]
    #                         eom = ep[om]
    #                         om_nutrient_index = list(eom).index(1.)
    #                         om_nutrient = conc[tracer_map[om][om_nutrient_index]]
    #                     else:
    #                         phy = p[0]
    #                         ephy = ep[phy]
    #                         phyto_nutrient_index = list(ephy).index(1.)
    #                         phyto_nutrient = conc[tracer_map[abbrev][phyto_nutrient_index]]
                        
    #                     # Get concentration of nutrient
    #                     elem_c = ec[cons][0]
    #                     ind_c = ic[cons][0]
    #                     nutrient = conc[tracer_map[cons][ind_c]]

    #                     # Calculate maximum nutrient uptake
    #                     specific_affinity = ids.index("specific_affinity")
    #                     max_uptake = params[specific_affinity] * nutrient * phyto

    #                     # Extract nutrient quota
    #                     quota_index = cell_quota_ids.index("si")
    #                     # nutrient_quota = cell_quota_opt[quota_index]
    #                     nutrient_quota = cell_quota_max[quota_index]

    #                     # Intracellular missing amount of nutrient
    #                     max_photo_rate = ids.index("max_photo_rate")
    #                     luxury_storage = ids.index("luxury_storage")
    #                     missing = params[max_photo_rate] * temp_regulation_factor * ( params[luxury_storage] * nutrient_quota * phyto - phyto_nutrient )

    #                     # Nutrient uptake based on net assimilation of C
    #                     assim_uptake = params[luxury_storage] * nutrient_quota * net_primary_production

    #                     # Actual uptake of nutrient
    #                     actual_uptake = np.minimum(max_uptake, missing + assim_uptake)

    #                     upt_switch = switch(actual_uptake)

    #                     uptake_to_phyto = upt_switch * actual_uptake
    #                     uptake_to_om = -actual_uptake * (1. - upt_switch)

    #                     uptake = uptake_to_phyto

    #                 # Update d_dt
    #                 d_dt[tracer_map[cons][0]] -= elem_c * uptake
    #                 if len(p) > 1:  
    #                     for i in range(len(ephy)):  d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)
    #                     for j in range(len(eom)):   d_dt[tracer_map[om][j]] += eom[j] *  uptake_to_om
                            
    #                 else:
    #                     for i in range(len(ephy)):  d_dt[tracer_map[abbrev][i]] += ephy[i] *  (uptake - uptake_to_om)

    #     return uptake

