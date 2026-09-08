import os
import sys
import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, monod, nutrient_limitation, tracer_elements, temperature_dependence, switch
from fractions import Fraction
np.set_printoptions(precision=20)
class Bacteria():
    """
    
    """


    def __init__(self, abbrev, base_element, physical, reactions, **tracer):
        
        # Variales that will be used later ---------------------------------------------------------------
        num_layers = physical["water_column"]["num_layers"]
        iters = physical["simulation"]["iters"]
        composition = physical["initial_conditions"][abbrev]["composition"]     # Initial concentrations
        if "scale" in physical["initial_conditions"][abbrev]:   scale = physical["initial_conditions"][abbrev]["scale"]     # Scaling factor (if initial concentration is split between multiple bacterioplankton groups)
        else:   scale = 1.  # No scaling

        # Add important keys ---------------------------------------------------------------
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Oxygen inhibition
        if "oxygen_inhibition" in tracer["parameters"]:
            # Initialize oxygen inhibition factor to array of 1. (if bacterioplankton is temperature limited this will be updated later otherwise will stay as 1.)
            if num_layers > 1:  self.oxy_limitation_factor = np.ones(num_layers-1, dtype=np.float64)
            else:   self.oxy_limitation_factor = np.float64(1.)

            if "oxygen_limited" in tracer["parameters"]["oxygen_inhibition"]:   self.oxygen_limited = tracer["parameters"]["oxygen_inhibition"]["oxygen_limited"]
            else:   self.oxygen_limited = False

            if self.oxygen_limited:
                # Default Hill exponent to 1 if not included in parameter list
                if "exponent" not in tracer["parameters"]["oxygen_inhibition"]:     tracer["parameters"]["oxygen_inhibition"]["exponent"] = np.float64(1.)

                self.oxy_inhib_ids = List.empty_list(unicode_type)
                self.oxy_inhib_params = List.empty_list(float64[:])
                
                for key,val in tracer["parameters"]["oxygen_inhibition"].items():
                    self.oxy_inhib_ids.append(key)
                    if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    self.oxy_inhib_params.append(val)

        # Organic Matter Partition
        if "om_partition" in tracer["parameters"]:
            # Initialize dict type
            dict_type = types.DictType(types.unicode_type, types.float64)
            lst = List.empty_list(dict_type)
            self.om_partition = Dict.empty(key_type=unicode_type,value_type=dict_type)

            # Convert parameter dictionary to numba typedDicts
            for outer_key in tracer["parameters"]["om_partition"]:
                d = Dict.empty(key_type=types.unicode_type, value_type = types.float64)
                inner_dict = tracer["parameters"]["om_partition"][outer_key]

                for inner_key in inner_dict:
                    d[inner_key] = np.float64(inner_dict[inner_key])

                self.om_partition[outer_key] = d

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

                self.temp_reg_ids = List.empty_list(unicode_type)
                self.temp_reg_params = List.empty_list(float64[:])
                
                for key,val in tracer["parameters"]["temperature_regulation"].items():
                    self.temp_reg_ids.append(key)
                    if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    self.temp_reg_params.append(val)
        else: self.temp_limited = False

        # Excretion
        if "excretion" in tracer["parameters"]:
            excr_dict_type = types.DictType(types.unicode_type, types.float64)
            self.excretion_ids = List.empty_list(unicode_type)  # stores produced tracer names (excretion parameters will be saved for each tracer individually)
            self.excretion_params = List.empty_list(excr_dict_type)

            for outer_key,inner_dict in tracer["parameters"]["excretion"].items():   # parse dictionary (key = tracer, val = inner dictionary)
                # Create temporary dictionary
                temp = Dict.empty(key_type=unicode_type, value_type=float64)

                # Add inner values to temporary dicitonary
                for inner_key,inner_val in inner_dict.items():
                    # Create numeric codes for excretion option string
                    if inner_key == "option":
                        if inner_val == "activity":     temp["option"] = 1.
                        elif inner_val == "constant":   temp["option"] = 2.
                        elif inner_val == "relaxation": temp["option"] = 3.

                    # Other inner values are already ints/floats
                    else:
                        temp[inner_key] = np.float64(inner_val) # conversion to make sure all values are floats

                # Add outer_key to excretion_ids and temp to excretion_params
                self.excretion_ids.append(outer_key)
                self.excretion_params.append(temp)

        # Mortality
        if "mortality" in tracer["parameters"]:
            # Create float option numbers for use in numba typed.List
            if "oxygen_limited" in tracer["parameters"]["mortality"]:
                # [0] False, [1] True
                if tracer["parameters"]["mortality"]["oxygen_limited"]:     tracer["parameters"]["mortality"]["oxygen_limited"] = 1
                else:   tracer["parameters"]["mortality"]["oxygen_limited"] = 0
            else:   tracer["parameters"]["mortality"]["oxygen_limited"] = 0     # Default to False if not in parameter list

            # Create list of mortality rates
            mort_rate = []  # [linear,quadratic]
            if "linear" in tracer["parameters"]["mortality"]["mortality_rate"]: mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["linear"]))
            else:   mort_rate.append(np.float64(0.))
            if "quadratic" in tracer["parameters"]["mortality"]["mortality_rate"]: mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["quadratic"]))
            else:   mort_rate.append(np.float64(0.))

            tracer["parameters"]["mortality"]["mortality_rate"] = np.array(mort_rate,dtype=np.float64)

            # Create list for temperature limitation
            temp_lim = []   # [linear,quadratic,oxygen]
            if "temp_limitation" in tracer["parameters"]["mortality"]:
                # [0] False, [1] True
                if "linear" in tracer["parameters"]["mortality"]["temp_limitation"]:    
                    if tracer["parameters"]["mortality"]["temp_limitation"]["linear"] == True:  temp_lim.append(1)
                    else:   temp_lim.append(0)
                else:   temp_lim.append(0)  # Default to False
                if "quadratic" in tracer["parameters"]["mortality"]["temp_limitation"]:    
                    if tracer["parameters"]["mortality"]["temp_limitation"]["quadratic"] == True:  temp_lim.append(1)
                    else:   temp_lim.append(0)
                else:   temp_lim.append(0)  # Default to False
                if "oxygen" in tracer["parameters"]["mortality"]["temp_limitation"]:    
                    if tracer["parameters"]["mortality"]["temp_limitation"]["oxygen"] == True:  temp_lim.append(1)
                    else:   temp_lim.append(0)
                else:   temp_lim.append(0)  # Default to False
            else:   temp_lim = [0,0,0]  # Set all to False

            tracer["parameters"]["mortality"]["temp_limitation"] = np.array(temp_lim,dtype=np.float64)

            self.mortality_ids = List.empty_list(unicode_type)
            self.mortality_params = List.empty_list(float64[:])
            
            for key,val in tracer["parameters"]["mortality"].items():
                self.mortality_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.mortality_params.append(val)

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
            if "convert_hs" in tracer["parameters"]["respiration"]:
                if isinstance(tracer["parameters"]["respiration"]["convert_hs"],str):
                    tracer["parameters"]["respiration"]["convert_hs"] = np.array([Fraction(tracer["parameters"]["respiration"]["convert_hs"])],dtype=np.float64)
                else:
                    tracer["parameters"]["respiration"]["convert_hs"] = np.array([tracer["parameters"]["respiration"]["convert_hs"]],dtype=np.float64)
            
            self.respiration_ids = List.empty_list(unicode_type)
            self.respiration_params = List.empty_list(float64[:])

            for key,val in tracer["parameters"]["respiration"].items():
                self.respiration_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.respiration_params.append(val)

        # Uptake
        if "uptake" in tracer["parameters"]:
            # Store uptake option string
            self.uptake_option = tracer["parameters"]["uptake"]["option"]

            # Create typed.Dict of uptake parameters and ids (will remain empty if option == "potential")
            self.uptake_ids = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(unicode_type))
            self.uptake_params = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(unicode_type))

            # Create typed.Dict of substrates (will remain empty if option != "balanced_substrate")
            self.substrates = Dict.empty(key_type=types.unicode_type, value_type=types.DictType(unicode_type,float64))
            self.max_growth_rate = 0.   # may be updated later
            self.half_sat_uptake = 0.   # may be updated later

            # Nested typed.Dict for coupled uptake (will remain empty if option == "potential")
            self.coupled_uptake = Dict.empty(
                key_type=types.unicode_type, 
                value_type=types.DictType(types.unicode_type, types.ListType(types.unicode_type))
            )
            
            # Create typed.Dict of uptake potentials (will remain empty if option == "direct")
            self.uptake_potential_rich = Dict.empty(key_type=types.unicode_type, value_type=float64)
            self.uptake_potential_poor = Dict.empty(key_type=types.unicode_type, value_type=float64)

            # Boolean to determine substrate correction (used if option == "potential")
            if "substrate_correction" in tracer["parameters"]["uptake"]:    self.substrate_correction = tracer["parameters"]["uptake"]["substrate_correction"]
            else:   self.substrate_correction = False

            if self.uptake_option == "balanced_substrate":

                self.max_growth_rate = tracer["parameters"]["uptake"]["max_growth_rate"]
                # self.half_sat_uptake = tracer["parameters"]["uptake"]["half_sat_uptake"]

                # for name,values in tracer["parameters"]["uptake"]["substrates"].items():
                #     temp = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

                #     for inner_key,inner_val in values.items():
                #         # Create numeric key for numerator
                #         if inner_key == "numerator":
                #             if inner_val == "self":         inner_val = 1.  # numerator is the consumed tracer
                #             elif inner_val == "substrate":  inner_val = 2.  # numerator is the substrate

                #     # Add temporary dictionary to substrates dictionary
                #     self.substrates[name] = temp


            if self.uptake_option == "potential":
                for key,val in tracer["parameters"]["uptake"]["potential_rich"].items():    self.uptake_potential_rich[key] = np.float64(val)
                for key,val in tracer["parameters"]["uptake"]["potential_poor"].items():    self.uptake_potential_poor[key] = np.float64(val)
                self.max_growth_rate = tracer["parameters"]["uptake"]["max_growth_rate"]

            elif self.uptake_option == "direct":
                # Parse keys in uptake parameters dictionary
                for outer_key,inner_dict in tracer["parameters"]["uptake"].items():
                    # Create temporary dictionary
                    # temp = Dict.empty(key_type=unicode_type, value_type=float64)
                    temp = Dict.empty(key_type=unicode_type, value_type=unicode_type)

                    # Skip if not a tracer
                    if outer_key in {"option", "substrate_correction", "potential_rich", "potential_poor"}:   continue

                    # Add uptake parameters for tracers to the uptake dictionaries
                    if "strategy" in inner_dict:
                        if inner_dict["strategy"] == "independent":
                            temp["strategy"] = inner_dict["strategy"]   # already a str, don't convert
                            # Convert int/float/bool to str
                            temp["max_growth_rate"] = str(inner_dict["max_growth_rate"])
                            temp["basal_metabolic_rate"] = str(inner_dict["basal_metabolic_rate"])
                            temp["max_growth_efficiency"] = str(inner_dict["max_growth_efficiency"])
                            temp["half_sat"] = str(inner_dict["half_sat"])
                            if "convert_uptake" in inner_dict:  
                                # Translate from fraction string to float64 (if necessary)
                                if isinstance(inner_dict["coupled_uptake"]["convert_uptake"],str):
                                    frac = Fraction(inner_dict["coupled_uptake"]["convert_uptake"])
                                    frac = np.float64(frac)
                                temp["convert_uptake"] = str(frac)

                        elif inner_dict["strategy"] == "coupled":
                            temp["strategy"] = inner_dict["strategy"]

                            # Create temporary dictionary of linked uptake tracers
                            temp_coupled_uptake = Dict.empty(key_type=unicode_type, value_type=types.ListType(unicode_type))

                            # Create list of coupled uptake links
                            links_list = List.empty_list(unicode_type)
                            for linked_tracer in inner_dict["link"]:    links_list.append(linked_tracer)
                            temp_coupled_uptake["links"] = links_list

                            # Empty list for method of linked calculation (max,min,sum,product)
                            method_list = List.empty_list(unicode_type)
                            if "method" in inner_dict["coupled_uptake"]:    method_list.append(inner_dict["coupled_uptake"]["method"])  # use provided method
                            else:   method_list.append("min")   # use min if method not provided
                            temp_coupled_uptake["method"] = method_list

                            # Translate from fraction string to float64 (if necessary)
                            if isinstance(inner_dict["coupled_uptake"]["convert_uptake"],str):
                                frac = Fraction(inner_dict["coupled_uptake"]["convert_uptake"])
                                frac = np.float64(frac)
                                # Convert decimal to str for list
                                temp_coupled_uptake["convert_uptake"] = List(str(frac))
                            else:
                                temp_coupled_uptake["convert_uptake"] = List(str(frac))

                            # Add to overall coupled uptake typed.Dict
                            self.coupled_uptake[outer_key] = temp_coupled_uptake

                        # Create temporary lists for ids and parameters
                        temp_ids = List.empty_list(unicode_type)
                        temp_params = List.empty_list(unicode_type)
                        for temp_key,temp_val in temp.items():
                            temp_ids.append(temp_key)
                            temp_params.append(temp_val)

                        # Add ids and parameters to uptake dictionaries
                        self.uptake_ids[outer_key] = temp_ids
                        self.uptake_params[outer_key] = temp_params


            # It might be simpler to NOT use numba typed dicts for uptake calculations to avoid constant conversions back and forward between string and floats but idk we'll see
            # If I were to use typed dicts, would need to convert "max_growth_rate", "basal_metabolic_rate", "max_growth_efficiency", and "half_sat" from floats to strings here
            # and from strings to floats in the uptake function

            # Excretion
        
        # Substrates
        if self.uptake_option == "balanced_substrate":
            self.substrate_ids = []
            self.substrate_params = []

            for element,substrate_info in tracer["parameters"]["substrates"].items():
                self.substrate_ids.append(element)
                self.substrate_params.append(substrate_info)
        
        # Uptake/Release
        if "uptake_release" in tracer["parameters"]:
            excr_dict_type = types.DictType(types.unicode_type, types.float64)
            self.uptake_release_ids = List.empty_list(unicode_type)  # stores produced tracer names (excretion parameters will be saved for each tracer individually)
            self.uptake_release_params = List.empty_list(excr_dict_type)

            for outer_key,inner_dict in tracer["parameters"]["uptake_release"].items():   # parse dictionary (key = tracer, val = inner dictionary)
                # Create temporary dictionary
                temp = Dict.empty(key_type=unicode_type, value_type=float64)

                # Add inner values to temporary dicitonary
                for inner_key,inner_val in inner_dict.items():
                    # Create numeric codes for excretion option string
                    temp[inner_key] = np.float64(inner_val) # conversion to make sure all values are floats

                # Add outer_key to excretion_ids and temp to excretion_params
                self.uptake_release_ids.append(outer_key)
                self.uptake_release_params.append(temp)

        # Add concentrations ---------------------------------------------------------------
        self.composition = List.empty_list(unicode_type)
        conc = []
        if len(composition) < 1:
            sys.exit("Detritus: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            # Reorder "composition" so "base_element" is at the start of the list
            for key in list(composition):
                if key != base_element: composition[key] = composition.pop(key)

            for key in list(composition):
                available_elements = ['c','n','p','fe']
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
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        # if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
        #     self.conc = np.zeros((len(self.composition),num_layers-1,iters),dtype=np.float64)
        #     for const in range(0,len(self.composition)):
        #         self.conc[const,:,0] = scale * conc[const][:-1] # Apply scaling factor here to prevent from applying multiple times in the above step
        # else:   # Model as single box
        #     self.conc = np.zeros((len(self.composition),num_layers,iters),dtype=np.float64)
        #     for const in range(0,len(self.composition)):
        #         self.conc[const,:,0] = scale * conc[const]      # Apply scaling factor here to prevent from applying multiple times in the above step
        # self.d_dt = np.zeros_like(self.conc[...,0],dtype=np.float64)
        # self.conc_ratio = np.ones_like(self.conc[...,0],dtype=np.float64)

        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.initial_conc = np.zeros((len(self.composition),num_layers-1),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.initial_conc[const,:] = scale * conc[const][:-1] # Apply scaling factor here to prevent from applying multiple times in the above step
        else:   # Model as single box
            self.initial_conc = np.zeros((len(self.composition),num_layers),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.initial_conc[const,:] = scale * conc[const]      # Apply scaling factor here to prevent from applying multiple times in the above step
        # self.d_dt = np.zeros_like(self.initial_conc[...],dtype=np.float64)
        # self.conc_ratio = np.ones_like(self.initial_conc[...],dtype=np.float64)

        # Add cell quotas ---------------------------------------------------------------
        # Create list of cell quota ids
        self.cell_quota_ids = List.empty_list(unicode_type)
        for element in self.composition:    self.cell_quota_ids.append(element)
        
        # for element in self.composition:    # Base element not in cell quotas (cell quota of base element would be 1.)
        #     if element != base_element: self.cell_quota_ids.append(element)
        
        # # Create list of maximum cell quotas
        # if "max" in tracer["parameters"]["cell_quota"]: 
        #     self.cell_quota_max = List.empty_list(float64)
        #     for element in self.cell_quota_ids: # Add quotas in same order as ids
        #         self.cell_quota_max.append(tracer["parameters"]["cell_quota"]["max"][element])

        # # Create list of minimum cell quotas
        # if "min" in tracer["parameters"]["cell_quota"]: 
        #     self.cell_quota_min = List.empty_list(float64)
        #     for element in self.cell_quota_ids: # Add quotas in same order as ids
        #         self.cell_quota_min.append(tracer["parameters"]["cell_quota"]["min"][element])

        # # Create list of optimal cell quotas
        # if "opt" in tracer["parameters"]["cell_quota"]: 
        #     self.cell_quota_opt = List.empty_list(float64)
        #     for element in self.cell_quota_ids: # Add quotas in same order as ids
        #         self.cell_quota_opt.append(tracer["parameters"]["cell_quota"]["opt"][element])
        
        if "cell_quota" in tracer["parameters"]:
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
        else:
            self.cell_quota_max = List.empty_list(float64)
            self.cell_quota_min = List.empty_list(float64)
            self.cell_quota_opt = List.empty_list(float64)


        # Add production arrays ---------------------------------------------------------------
        self.upt = Dict.empty(
                key_type=types.unicode_type, 
                value_type=types.float64[:]
            )   # Uptake

        # Add relevant reactions
        self.reactions = []
        for reac in reactions:
            # Add reaction to dictionary
            if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
            else:   consumed = {"empty": "empty"}
            if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
            else:   produced = {"empty": "empty"}
            if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
                self.reactions.append(reac)

            # Delete "loss" reactions if this tracer is produced
            if ( reac["type"] == "loss" ) and ( abbrev in produced.keys() ):    self.reactions.pop()

        # Reorder uptake reactions in case of coupled uptake
        for i in range(len(self.reactions)):
            if self.reactions[i]["type"] == "uptake":
                # Get element of nutrient being consumed
                nutrient_element = list(self.reactions[i]["consumed"].values())[0][0]

                # Bump reaction to end of list if the nutrient is coupled to a different uptake rate
                if nutrient_element in self.coupled_uptake: self.reactions.append(self.reactions.pop(i))
        
        # Reorder reactions (uptake needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]

        # Boolean to determine whether activity and basal respiration will be calculated
        self.calc_respiration = False
        for reac in self.reactions:
            if reac["type"] == "respiration":
                self.calc_respiration = True
                break


    # def bac(self, iter, base_element, physical, tracers):
    def bac(self, base_element, temperature, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers):

        # Zero out variables
        for nut in self.upt:                                # Uptake
            self.upt[nut] = np.zeros_like(conc[0],dtype=np.float64)
        
        # Calculate oxygen limitation factor (if necessary)
        if self.oxygen_limited:
            # # Optional minimum concentration for aerobic/anerobic operations
            # if "min_o2" in self.oxy_inhib_params:  o2 = np.maximum(tracers["o2"].conc[...,iter] , self.oxygen_inhibition["min_o2"] * np.ones_like(tracers["o2"].conc[...,iter]))
            # else:   o2 = tracers["o2"].conc[...,iter]
            # Optional minimum concentration for aerobic/anerobic operations
            if "min_o2" in self.oxy_inhib_ids:  o2 = np.maximum(conc[tracer_map["o2"][0]] , self.oxy_inhib_params[self.oxy_inhib_ids.index("min_o2")] * np.ones_like(conc[tracer_map["o2"][0]]))
            else:   o2 = conc[tracer_map["o2"][0]]
            self.oxy_limitation_factor = monod(o2, self.oxy_inhib_params[self.oxy_inhib_ids.index("half_sat")], self.oxy_inhib_params[self.oxy_inhib_ids.index("exponent")])

        # Calculate temp regulation factor (if necessary)
        if self.temp_limited:
            self.temp_regulation_factor = temperature_dependence(temperature, self.temp_reg_ids, self.temp_reg_params)

        # Calculate potential uptake (if necessary)
        if self.uptake_option == "potential":
            actual_uptake, realized_uptake, base_uptake = self.calculate_realized_uptake(base_element, self.max_growth_rate, conc, conc_ratio, tracer_map, tracers)
        else:
            actual_uptake = np.zeros_like(conc[tracer_map[self.abbrev][0]])
            realized_uptake = np.zeros_like(conc[tracer_map[self.abbrev][0]])
            base_uptake = self.upt.copy()

        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "excretion":     
                if self.uptake_option == "direct":  # "actual_uptake" is sum off all uptakes in self.upt if using "direct" option, otherwise use calculated "actual_uptake" from above
                    for uptake_array in self.upt.values():  actual_uptake += uptake_array
                self.excretion(base_element, c, p, ec, ep, ic, ip, self.cell_quota_ids, self.cell_quota_opt, self.excretion_ids, self.excretion_params, actual_uptake, conc, conc_ratio, d_dt, tracer_map, tracer_type, self.composition)
            if reac["type"] == "mortality":     
                om_composition = tracers[p[0]].composition
                bac_composition = self.composition
                self.mortality(c, p, ec, ep, ic, ip, self.mortality_ids, self.mortality_params, self.om_partition, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, bac_composition, om_composition)
            if reac["type"] == "respiration":
                bact_limitation_factor = self.respiration(self.abbrev, base_element, c, p, ec, ep, ic, ip, self.respiration_ids, self.respiration_params, actual_uptake, self.oxy_limitation_factor, self.temp_regulation_factor, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "uptake":
                if self.uptake_option == "balanced_substrate":
                    self.calculate_balanced_substrate_uptake(c, ec, p, ep, self.upt, conc, d_dt, tracer_map, tracers)
                else:
                    bac_composition = self.composition
                    om_composition = tracers[c[0]].composition
                    self.uptake(self.abbrev, c, p, ec, ep, ic, ip, self.uptake_ids, self.uptake_params, self.uptake_option, self.upt, actual_uptake, realized_uptake, base_uptake, self.coupled_uptake, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, bac_composition, om_composition)
            if reac["type"] == "uptake_release":
                self.uptake_release(base_element, c, p, ec, ep, ic, ip, self.cell_quota_ids, self.cell_quota_opt, self.uptake_release_ids, self.uptake_release_params, conc, conc_ratio, d_dt, tracer_map, tracer_type, self.composition)
        
        # Set bacteria limitation factor to zero if no respiration
        if not self.calc_respiration:   
            bact_limitation_factor = 0.

        return bact_limitation_factor


    def add_nutrient(self, nutrients):
        """
        Add "nutrients" to bacterioplankton and append dictionary of uptake rates. 
        "Nutrients" are organic matter pools, only one array is used to store uptake of base element from each organic matter pool.
        """
        # zeros = List.empty_list(float64[:])
        # zeros.append(np.zeros(self.conc.shape[1],dtype=np.float64))
        # for nut in nutrients:
        #     self.upt[nut] = np.zeros_like(self.conc[0,:,0],dtype=np.float64)
        zeros = List.empty_list(float64[:])
        zeros.append(np.zeros(self.initial_conc.shape[1],dtype=np.float64))
        for nut in nutrients:
            self.upt[nut] = np.zeros_like(self.initial_conc[0,:],dtype=np.float64)


    def calculate_balanced_substrate_uptake(self, c, ec, p, ep, upt, conc, d_dt, tracer_map, tracers):

        cons = c[0]
        elem_c = list(ec[cons])
        element = tracers[cons].composition[elem_c.index(1.)]

        prod = p[0]
        elem_p = ep[prod]

        bac_element_index = self.substrate_ids.index(element)

        # Extract substrate information for the element being consumed
        element_substrate_dict = self.substrate_params[bac_element_index]

        components = np.zeros((len(self.substrate_params[bac_element_index])-2,len(self.initial_conc[0])),dtype=np.float64)
        i=0
        for key,val in self.substrate_params[bac_element_index].items():
            if key not in {"denominator","half_sat_uptake"}:
                index = tracers[key].composition.index(self.substrate_ids[bac_element_index])
                components[i] = conc[tracer_map[key][index]] * element_substrate_dict[key]["coefficient"]
                i += 1
            
            substrate = np.min(components,axis=0)

        if element_substrate_dict[cons]["numerator"] == "self":           numerator = conc[tracer_map[cons][elem_c.index(1.)]]
        elif element_substrate_dict[cons]["numerator"] == "substrate":    numerator = substrate

        denominator = np.zeros_like(numerator)
        for key,val in element_substrate_dict["denominator"].items():
            if key == "substrate":  denominator += substrate
            else:   denominator += conc[tracer_map[key][tracers[key].composition.index(val[0])]]

        denominator += element_substrate_dict["half_sat_uptake"]

        uptake = self.max_growth_rate * (numerator / denominator) * conc[tracer_map[prod][self.composition.index(element)]]

        # Update d_dt
        for i in range(len(elem_c)):
            d_dt[tracer_map[cons][i]] -= elem_c[i] * uptake
        for j in range(len(elem_p)):
            d_dt[tracer_map[prod][j]] += elem_p[j] * uptake

        upt[cons] = uptake

    
    def calculate_realized_uptake(self, base_element, growth_rate, conc, conc_ratio, tracer_map, tracers):
        # Get concentration of base element in bacterioplankton
        base_index = self.composition.index(base_element)
        bac = conc[tracer_map[self.abbrev][base_index]]

        # List of organic matter pools for uptake
        om_pools = list(self.uptake_potential_rich.keys())

        if self.substrate_correction:   # Apply substrate correction
            # Calculate potential uptake rate
            potential_uptake = growth_rate * self.temp_regulation_factor * self.nutrient_colimitation_factor * bac

            # Initialize correction matrix
            correction = np.ones((len(om_pools),len(bac)),dtype=np.float64)

            # Calculate correction of substrate quality depending on nutrient content for each organic matter pool
            for om in om_pools:
                # Initialize correction to 1.
                om_correction = np.ones_like(bac)

                # Get composition of organic matter pool
                om_composition = tracers[om].composition

                # Use optimal bacteria cell quota to calculate correction
                for const in om_composition:
                    if const != base_element:   # ignore base element
                        const_quota = self.cell_quota_opt[const]    # cell quota of element in bacteria

                        om_index = tracers[om].composition.index(const)     # index of element in organic matter pool
                        om_ratio = conc_ratio[tracer_map[om][om_index]]     # concentration ratio of element in organic matter pool

                        om_correction = np.minimum(om_correction, (om_ratio/const_quota))

                # Add correction to matrix of correction values
                correction[om_pools.index(om)] = om_correction
            
        else:   # no substrate correction
            # Calculate potential uptake rate
            potential_uptake = growth_rate * self.temp_regulation_factor * bac

            # Set correction to 1.
            correction = np.ones((len(om_pools),len(bac)),dtype=np.float64)


        # Calcualate realized uptake
        realized_uptake = 1.E-20 * np.ones_like(bac)    # Initialize to 1.E-20 to prevent divide by 0.
        
        # base_uptake = {}    # Create uptake dictionary
        base_uptake = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])  # Create typed.Dict for uptake of base element in nutrient

        # Calculate potential uptake of base element for organic matter pools
        for om in om_pools:
            # Get composition of organic matter pool
            om_composition = tracers[om].composition

            # Get index of base element in organic matter pool
            om_base_index = tracers[om].composition.index(base_element)

            # Calculate substrate uptake of base element in each organic matter pool
            base_uptake[om] = ( ( self.uptake_potential_rich[om] * correction[om_pools.index(om)] ) + ( self.uptake_potential_poor[om] * ( 1. - correction[om_pools.index(om)] ) ) ) * conc[tracer_map[om][om_base_index]]
        
            # Add to total for realized uptake
            realized_uptake += base_uptake[om]

        # Calculate actual uptake (minimum of potential and realized uptakes)
        actual_uptake = np.minimum(potential_uptake, realized_uptake)
        
        return actual_uptake, realized_uptake, base_uptake


    @staticmethod
    @njit
    def excretion(base_element, c, p, ec, ep, ic, ip, cell_quota_ids, cell_quota_opt, excretion_ids, excretion_params, actual_uptake, conc, conc_ratio, d_dt, tracer_map, tracer_type, composition):
        """
        Definition:: Calculates excretion of bacterioplankton to nutrient pool.
                     Excretion can be represented either as a constant rate or as a fraction of bacterioplankton uptake on organic matter.
        Return:: Excretion rate
        """
        # Nutrient elements
        nutrient_elements = Dict.empty(key_type=types.unicode_type,value_type=types.unicode_type)
        nutrient_elements["no3"] = "n"
        nutrient_elements["nh4"] = "n"
        nutrient_elements["po4"] = "p"

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons]
        ind_p = ip[prod]
        
        # Get concentration of constituent in bacterioplankton
        # bac = conc[tracer_map[cons][0]]
        for element in range(len(elem_c)):
            if elem_c[element] == 1:
                idx = element
                break
        bac = conc[tracer_map[cons][list(elem_c).index(1)]]
        excreted_element = composition[list(elem_c).index(1)]
        base_index = composition.index(base_element)    # Index of base element

        # Extract parameter indices
        nutrient_excretion = excretion_ids.index(prod)

        # Determine type of "produced" nutrient (inorganic nutrient or organic matter pool)
        if tracer_type[tracer_map[prod][0]] == "inorganic": # Extract cell quota id of nutrient element if inorganic
            nutrient_index = composition.index(nutrient_elements[prod])     # Index of excreted nutrient
            quota_index = cell_quota_ids.index(nutrient_elements[prod])
            # external_nut_lim = nutrient_limitation(conc[tracer_map[prod][0]], excretion_params[nutrient_excretion]["half_sat"])
        # Uses minimum of all quotas if organic matter pool, calculated below

        # Excretion
        if excretion_params[nutrient_excretion]["option"] == 1.:    # based on activity
            excretion = actual_uptake * (1. - excretion_params[nutrient_excretion]["activity_respiration_frac"]) * (excretion_params[nutrient_excretion]["activity_respiration_frac"] * excretion_params[nutrient_excretion]["excretion_rate"])

        elif excretion_params[nutrient_excretion]["option"] == 2.:  # constant excretion rate
            excretion = excretion_params[nutrient_excretion]["excretion_rate"] * bac

        elif excretion_params[nutrient_excretion]["option"] == 3.:  # based on relaxation time
            # Calculate excess above optimal cell quota
            if tracer_type[tracer_map[prod][0]] == "inorganic":
                excess = conc_ratio[tracer_map[cons][nutrient_index]] - cell_quota_opt[quota_index]
            else:
                max_quota = np.zeros_like(bac)  # initialize array
                for const in composition:   # parse composition list
                    if const in cell_quota_ids:
                        const_index = composition.index(const)  # constituent index
                        quota_index = cell_quota_ids.index(const)   # cell quota index for constituent
                        hold = 1. - (conc_ratio[tracer_map[cons][const_index]]/cell_quota_opt[quota_index])
                        max_quota = np.maximum(max_quota, hold)

                excess = np.maximum(np.zeros_like(bac), max_quota)

            # Calculate excretion
            excretion = excess * excretion_params[nutrient_excretion]["relaxation_timescale"] * bac

        # Apply correction for excess excretion (if necessary)
        # if base_index == nutrient_index:    # not necessary
        # if base_element == excreted_element:    # not necessary
        #     pass
        # else:   # excess correction
        #     element_ratio = bac / conc[tracer_map[cons][base_index]]
        #     excretion = excretion * np.maximum(0., element_ratio - cell_quota_opt[quota_index])

        #     # Aply nutrient limitation factor to inorganic nutrient excretion
        #     if tracer_type[tracer_map[prod][0]] == "inorganic": excretion *= external_nut_lim


        # Update d_dt
        # d_dt[tracer_map[cons][nutrient_index]] -= excretion
        for i in range(len(tracer_map[cons])):  d_dt[tracer_map[cons][i]] -= elem_c[i] * excretion
        for j in range(len(tracer_map[prod])):  d_dt[tracer_map[prod][j]] += elem_p[j] * excretion

        return


    @staticmethod
    @njit
    def mortality(c, p, ec, ep, ic, ip, mortality_ids, mortality_params, om_partition, temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, bac_composition, om_composition):
        # d_0b = specific mortality rate (linear)
        # d_B_d = density specific mortality rate (quadratic)

        # Extract parameters
        mortality_rate = mortality_ids.index("mortality_rate")
        temp_limitation = mortality_ids.index("temp_limitation")

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons][0]
        ind_p = ip[prod][0]

        bac = conc[tracer_map[cons][ind_c]]

        # Calculate mortality rate (temperature regulated linear component)
        # mortality = ( mortality_params[mortality_rate][0] * temp_regulation_factor * bac ) + ( mortality_params[mortality_rate][1] * (bac**2) )

        linear = mortality_params[mortality_rate][0] * bac
        quadratic = mortality_params[mortality_rate][1] * (bac**2)

        # Apply temperature limitation if necessary (1. = True)
        if mortality_params[temp_limitation][0] == 1.:  linear *= temp_regulation_factor
        if mortality_params[temp_limitation][1] == 1.:  quadratic *= temp_regulation_factor

        mortality = linear + quadratic

        # Calculate concentration ratios
        ratios = np.zeros((len(bac_composition),len(bac)),dtype=np.float64)
        for const in bac_composition:
            if const in om_composition:
                index_bac = bac_composition.index(const)
                index_om = om_composition.index(const)
                ratios[index_om] = conc_ratio[tracer_map[cons][index_bac]]

        # Update d_dt
        if prod in om_partition:    # Apply partition (if necessary)
            for const in om_partition[prod]:
                om_const_index = om_composition.index(const)    # Get index of constituent in organic matter pool
                if const in bac_composition:    # Only apply rate if constituent is also in zooplankton
                    bac_const_index = bac_composition.index(const)  # Get index of constituent in zooplankton
                    d_dt[tracer_map[cons][bac_const_index]] -= elem_c[bac_const_index] * conc_ratio[tracer_map[cons][bac_const_index]] * mortality * om_partition[prod][const]
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * mortality * om_partition[prod][const]

        else:
            for const in om_composition:
                om_const_index = om_composition.index(const)
                if const in bac_composition:
                    bac_const_index = bac_composition.index(const)
                    d_dt[tracer_map[cons][bac_const_index]] -= elem_c[bac_const_index] * conc_ratio[tracer_map[cons][bac_const_index]] * mortality
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * mortality

        return
    
    @staticmethod
    @njit
    def respiration(abbrev, base_element, c, p, ec, ep, ic, ip, respiration_ids, respiration_params, actual_uptake, oxy_limitation_factor, temp_regulation_factor, conc, d_dt, tracer_map, composition):

        # Extract parameter indices
        activity_respiration_frac = respiration_ids.index("activity_respiration_frac")
        basal_respiration_rate = respiration_ids.index("basal_respiration_rate")
        anoxic_respiration_frac = respiration_ids.index("anoxic_respiration_frac")
        conv_o2 = False
        if "convert_o2" in respiration_ids:
            conv_o2 = True
            convert_o2 = respiration_ids.index("convert_o2")
        conv_co2 = False
        if "convert_co2" in respiration_ids:
            conv_co2 = True
            convert_co2 = respiration_ids.index("convert_co2")
        conv_hs = False
        if "convert_hs" in respiration_ids:
            conv_hs = True
            convert_hs = respiration_ids.index("convert_hs")

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

        # Get base element concentration
        bac = conc[tracer_map[abbrev][base_index]]

        # Respiration
        activity_respiration = ( respiration_params[activity_respiration_frac] + respiration_params[anoxic_respiration_frac] * (1. - oxy_limitation_factor) ) * actual_uptake
        basal_respiration = respiration_params[basal_respiration_rate] * temp_regulation_factor * bac
        respiration = activity_respiration + basal_respiration

        # Update d_dt
        d_dt[tracer_map[abbrev][base_index]] -= respiration
        if "o2" in c:   
            if conv_o2:     d_dt[tracer_map["o2"][0]] -= oxy_limitation_factor * respiration * respiration_params[convert_o2]
            else:           d_dt[tracer_map["o2"][0]] -= oxy_limitation_factor * respiration

        if p:   # "produced" not an empty list
            for element in range(len(p)):
                if p[element] == "co2":
                    if base_element == "c":     d_dt[tracer_map["co2"][0]] += respiration  
                    else:                       
                        if conv_co2:    d_dt[tracer_map["co2"][0]] += respiration * respiration_params[convert_co2]
                        else:           d_dt[tracer_map["co2"][0]] += respiration
                elif p[element] == "hs":
                    if conv_hs:         d_dt[tracer_map["hs"][0]] += respiration * (1. - oxy_limitation_factor) * respiration_params[convert_hs]
                    else:               d_dt[tracer_map["hs"][0]] += respiration * (1. - oxy_limitation_factor)
                else:   
                    d_dt[tracer_map[p[element]][ip[p[element]][0]]] += respiration

        # Calculate bacteria limiting factor for denitrification (if necessary)
        if conv_hs: bact_limitation_factor = (1. - oxy_limitation_factor) * respiration * respiration_params[convert_hs]
        else:       bact_limitation_factor = (1. - oxy_limitation_factor) * respiration 

        return bact_limitation_factor


    @staticmethod
    # @njit
    def uptake(abbrev, c, p, ec, ep, ic, ip, uptake_ids, uptake_params, uptake_option, upt, actual_uptake, realized_uptake, base_uptake, coupled_uptake_dict, temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, bac_composition, om_composition):
        
        # Extract dict
        cons = c[0] # organic matter
        prod = p[0] # bacterioplankton
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons][0]
        ind_p = ip[prod][0]

        # Extract concentration ratio of consumed "nutrient"
        ratios = np.zeros((len(bac_composition),len(conc[tracer_map[abbrev][0]])))
        for const in bac_composition:
            if const in om_composition:
                index_om = om_composition.index(const)
                index_bac = bac_composition.index(const)
                ratios[index_bac] = conc_ratio[tracer_map[cons][index_om]]
        
        if uptake_option == "direct":

            ids = uptake_ids[cons]
            params = uptake_params[cons]
        
            # Identify strategy for uptake rate calculation
            strategy = ids.index("strategy")

            if params[strategy] == "independent":    # "independent" uptake
                # Extract parameter ids
                max_growth_rate_id = ids.index("max_growth_rate")
                basal_metabolic_rate_id = ids.index("basal_metabolic_rate")
                max_growth_efficiency_id = ids.index("max_growth_efficiency")
                half_sat_id = ids.index("half_sat")

                # Parameters stored as strings, convert to floats
                max_growth_rate = np.float64(params[max_growth_rate_id])
                basal_metabolic_rate = np.float64(params[basal_metabolic_rate_id])
                max_growth_efficiency = np.float64(params[max_growth_efficiency_id])
                half_sat = np.float64(params[half_sat_id])

                convert = False
                if "convert_uptake" in ids:  
                    convert = True
                    convert_uptake = np.float64(params[ids.index("convert_uptake")])

                # Get concentration of base element
                bac = conc[tracer_map[abbrev][ind_p]]
                om = conc[tracer_map[cons][ind_c]]

                # Calculate maximum uptake rate
                # max_uptake = ( params[max_growth_rate] + params[basal_metabolic_rate] ) / params[max_growth_efficiency]
                max_uptake = ( max_growth_rate + basal_metabolic_rate ) / max_growth_efficiency

                # Calculate actual uptake
                # uptake = temp_regulation_factor * max_uptake * monod(om, params[half_sat], 1.) * bac
                uptake = temp_regulation_factor * max_uptake * monod(om, half_sat, 1.) * bac

                if convert: uptake *= convert_uptake

                # Update d_dt
                for i in range(len(elem_c)):    d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * uptake
                for j in range(len(elem_p)):    d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * uptake

                # Update uptake rate dictionary
                upt[cons] = uptake

            elif params[strategy] == "coupled":      # "coupled" uptake
                coupled_uptake = coupled_uptake_dict[cons]
                linked_nutrients = coupled_uptake["links"]
                convert = False
                if "convert_uptake" in coupled_uptake:  convert = True

                # Extract uptake rates of linked nutrients
                uptake_rates = List.empty_list(float64[:])

                for nut in linked_nutrients:
                    uptake_rates.append(upt[nut])    # uptake rate of base element in organic matter pool

                # Calculate total linked uptake rate if multiple linked nutrients are used
                linked_uptake = uptake_rates[0].copy()
                if linked_nutrients and len(linked_nutrients) > 1:   # use numpy "maximum" to ensure minimum uptake of 0.
                    if coupled_uptake["method"] == "max":      # "max"
                        for i in range(1,len(uptake_rates)):
                            linked_uptake = np.maximum(linked_uptake,uptake_rates[i])

                    elif coupled_uptake["method"] == "min":    # "min"
                        for i in range(1,len(uptake_rates)):
                            linked_uptake = np.minimum(linked_uptake, uptake_rates[i])
                        
                    elif coupled_uptake["method"] == "sum":    # "sum"
                        for i in range(1,len(uptake_rates)):
                            linked_uptake += uptake_rates[i]

                    elif coupled_uptake["method"] == "product":    # "product"
                        for i in range(1,len(uptake_rates)):
                            linked_uptake *= uptake_rates[i]
                    
                # minimum uptake 0.
                linked_uptake = np.maximum(linked_uptake, np.zeros_like(uptake_rates[0]))

                if convert:
                    # convert "convert_uptake" from string to float
                    convert_uptake = np.float64(coupled_uptake["convert_uptake"][0])
                    uptake = np.zeros(len(conc[tracer_map[abbrev][0]]),dtype=np.float64)

                    for depth in range(len(conc[tracer_map[abbrev][0]])):
                        uptake[depth] = linked_uptake[depth] * convert_uptake

                # Update d_dt (apply concentration ratios)
                for i in range(len(elem_c)):    d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * uptake
                for j in range(len(elem_p)):    d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * uptake

                # Update uptake rate dictionary
                upt[cons] = uptake


        elif uptake_option == "potential":
            # Calculate uptake rate
            uptake = actual_uptake * base_uptake[cons] / realized_uptake

            # Update d_dt
            for i in range(len(elem_c)):    d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * uptake
            for j in range(len(elem_p)):    d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * uptake
            
            # Update uptake rate didctionary
            upt[cons] = uptake


    @staticmethod
    @njit
    def uptake_release(base_element, c, p, ec, ep, ic, ip, cell_quota_ids, cell_quota_opt, uptake_release_ids, uptake_release_params, conc, conc_ratio, d_dt, tracer_map, tracer_type, composition):
        """
        Definition:: Calculate the uptake/release between bacterioplankton and the nutrient pool.
                     
        """
        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons]
        ind_p = ip[prod]

        # Extract nutrient for parameters
        nutrient = uptake_release_ids.index(prod)
        # half_sat = uptake_release_ids.index("half_sat")
        # relaxation_timescale = uptake_release_ids.index("relaxation_timescale")
        
        # Get concentration of base element in bacterioplankton
        base_index = composition.index(base_element)    # Index of base element
        bac = conc[tracer_map[cons][base_index]]

        # Get cell quotas
        element_index = list(elem_c).index(1)
        element = composition[element_index]
        current_ratio = conc_ratio[tracer_map[cons][element_index]]  # Current cell ratio
        optimal_quota = cell_quota_opt[cell_quota_ids.index(element)]   # Extract the optimal cell quota of the nutrient being release

        # Calculate external nutrient limitation
        nut_lim = nutrient_limitation(conc[tracer_map[prod][0]], uptake_release_params[nutrient]["half_sat"])
        
        # Calculate uptake/release
        upt_rel = (current_ratio - optimal_quota) * uptake_release_params[nutrient]["relaxation_timescale"] * bac

        # Calculate switch vectors for direction between bacteria and nutrient
        pos_switch = switch(upt_rel)
        neg_switch = switch(-upt_rel)

        # # Update d_dt
        # d_dt[tracer_map[cons][element_index]] -= upt_rel * (pos_switch - (nut_lim * neg_switch))    # To bacteria
        # d_dt[tracer_map[prod][0]] += upt_rel * (pos_switch - (nut_lim * neg_switch))                # To nutrient

        # Update d_dt
        d_dt[tracer_map[cons][element_index]] -= upt_rel * (pos_switch + (nut_lim * neg_switch))    # To bacteria
        d_dt[tracer_map[prod][0]] += upt_rel * (pos_switch + (nut_lim * neg_switch))                # To nutrient

        return
    