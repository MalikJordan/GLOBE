import copy
import os
import sys
import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from functions.other_functions import concentration_ratio, monod, nutrient_limitation, string_to_float, temperature_dependence, tracer_elements
from fractions import Fraction
np.set_printoptions(precision=20)
class Zooplankton():
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
                available_elements = ['c','n','p']
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
                                    conc.append( composition[key][factor_element] * conc[index] * np.ones(num_layers,dtype=np.float64)) 

                else:
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
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
        
        # Create list of optimal cell quotas
        if "opt" in tracer["parameters"]["cell_quota"]: 
            self.cell_quota_opt = List.empty_list(float64)
            for element in self.cell_quota_ids: # Add quotas in same order as ids
                self.cell_quota_opt.append(tracer["parameters"]["cell_quota"]["opt"][element])

        # Add efficiencies ---------------------------------------------------------------
        if isinstance(tracer["parameters"]["efficiency"]["assimilation"],str):
            frac = Fraction(tracer["parameters"]["efficiency"]["assimilation"])
            frac = np.float64(tracer["parameters"]["efficiency"]["assimilation"])
            self.assimilation_efficiency = frac
        else:   self.assimilation_efficiency = np.float64(tracer["parameters"]["efficiency"]["assimilation"])

        if isinstance(tracer["parameters"]["efficiency"]["ingestion"],str):
            frac = Fraction(tracer["parameters"]["efficiency"]["ingestion"])
            frac = np.float64(tracer["parameters"]["efficiency"]["ingestion"])
            self.ingestion_efficiency = frac
        else:   self.ingestion_efficiency = np.float64(tracer["parameters"]["efficiency"]["ingestion"])

        # Add rate parameters ---------------------------------------------------------------
        # Egestion
        # if "egestion" in tracer["parameters"]:
        #     for key in tracer["parameters"]["egestion"]:
        #         if isinstance(tracer["parameters"]["egestion"][key],str):   # calculated based on grazing rates
        #             tracer["parameters"]["egestion"][key] = -1
        #         else:   # constant excretion rate
        #             tracer["parameters"]["egestion"][key] = np.float64(tracer["parameters"]["egestion"][key])

        #     self.egestion_ids = List.empty_list(unicode_type)
        #     self.egestion_params = List.empty_list(float64)

        #     for key,val in tracer["parameters"]["egestion"].items():
        #         self.egestion_ids.append(key)
        #         if not isinstance(val, np.float64): val = np.float64(val)
        #         self.egestion_params.append(val)

        if "egestion" in tracer["parameters"]:
            self.unassimilated_egestion = List.empty_list(unicode_type)
            for const in tracer["parameters"]["egestion"]["unassimilated"]: self.unassimilated_egestion.append(const)
            # self.unassimilated_egestion = tracer["parameters"]["egestion"]["unassimilated"]
        else:
            self.unassimilated_egestion = List.empty_list(unicode_type)

        # Excretion
        if "excretion" in tracer["parameters"]:
            # Create destination and constituent lists
            constituents = List.empty_list(unicode_type)
            destination = List.empty_list(unicode_type)
            for key,val in tracer["parameters"]["excretion"]["destination"].items():
                constituents.append(key)
                destination.append(val)

            # Create function list
            function = List.empty_list(unicode_type)
            for const in constituents:
                # function can include numeric value for constant excretion rate, convert to string for typed.List (will be converted back inside of excretion calculation)
                if not isinstance(tracer["parameters"]["excretion"]["function"][const],str):    function.append(str(tracer["parameters"]["excretion"]["function"][const]))
                else:   function.append(tracer["parameters"]["excretion"]["function"][const])

            self.excretion_ids = List(["constituents","destination","function"])
            self.excretion_params = List.empty_list(types.ListType(unicode_type))
            self.excretion_params.append(constituents)
            self.excretion_params.append(destination)
            self.excretion_params.append(function)

            # for key in tracer["parameters"]["excretion"]:
            #     if isinstance(tracer["parameters"]["excretion"][key],str):   # calculated based on grazing rates
            #         tracer["parameters"]["excretion"][key] = -1
            #     else:   # constant excretion rate
            #         tracer["parameters"]["excretion"][key] = np.float64(tracer["parameters"]["excretion"][key])

            # self.excretion_ids = List.empty_list(unicode_type)
            # self.excretion_params = List.empty_list(float64)

            # for key,val in tracer["parameters"]["excretion"].items():
            #     self.excretion_ids.append(key)
            #     if not isinstance(val, np.float64): val = np.float64(val)
            #     self.excretion_params.append(val)

        # Grazing 
        if "grazing" in tracer["parameters"]:
            # Create float option numbers for use in numba typed.List
            if "function" in tracer["parameters"]["grazing"]:
                if tracer["parameters"]["grazing"]["function"] == "holling-1":      tracer["parameters"]["grazing"]["function"] = 1
                elif tracer["parameters"]["grazing"]["function"] == "holling-2":    tracer["parameters"]["grazing"]["function"] = 2
                elif tracer["parameters"]["grazing"]["function"] == "holling-3":    tracer["parameters"]["grazing"]["function"] = 3
                elif tracer["parameters"]["grazing"]["function"] == "ivlev":        tracer["parameters"]["grazing"]["function"] = 4

            # Create float option for use of capture efficiency in grazing calculations (default to "True")
            # Typically "True" for microzooplankton and "False" for mesozooplankton (grazing based on clearance rate rather than capture efficiency)
            if "use_capture_efficiency" in tracer["parameters"]["grazing"]:
                if tracer["parameters"]["grazing"]["use_capture_efficiency"] == True:   tracer["parameters"]["grazing"]["use_capture_efficiency"] = 1.
                else:   tracer["parameters"]["grazing"]["use_capture_efficiency"] = 0.
            else:   # Default to True
                tracer["parameters"]["grazing"]["use_capture_efficiency"] = 1.

            # Create float option for limiting factor (default to "half_saturation")
            # Typically "half_saturation" for microzooplankton and "clearance_rate" for mesozooplankton
            if "feeding_model" in tracer["parameters"]["grazing"]:
                if tracer["parameters"]["grazing"]["feeding_model"] == "clearance_rate":   tracer["parameters"]["grazing"]["feeding_model"] = 1
                elif tracer["parameters"]["grazing"]["feeding_model"] == "half_saturation":      tracer["parameters"]["grazing"]["feeding_model"] = 2
            else:   # Default to "half_saturation"
                tracer["parameters"]["grazing"]["feeding_model"] = 2

            self.grazing_ids = List.empty_list(unicode_type)
            self.grazing_params = List.empty_list(float64)

            for key,val in tracer["parameters"]["grazing"].items():
                if key == "grazing_preferences":    pass
                else:
                    self.grazing_ids.append(key)
                    self.grazing_params.append(np.float64(val))
               
            self.grazing_preferences = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
            for key,val in tracer["parameters"]["grazing"]["grazing_preferences"].items():
                self.grazing_preferences[key] = np.float64(tracer["parameters"]["grazing"]["grazing_preferences"][key])

            self.prey_availability = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
            self.grazing_rates = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])

        # Metabolic Release
        if "metabolic_release" in tracer["parameters"]:
            # set basal metabolic rate to 0. if not included in parameter list
            if "basal_metabolic_rate" not in tracer["parameters"]["metabolic_release"]: self.basal_metabolic_rate = 0.
            else: self.basal_metabolic_rate = tracer["parameters"]["metabolic_release"]["basal_metabolic_rate"]

            # Create destination and constituent lists
            constituents = List.empty_list(unicode_type)
            destination = List.empty_list(unicode_type)
            for key,val in tracer["parameters"]["metabolic_release"]["destination"].items():
                constituents.append(key)
                destination.append(val)

            # Create function list
            function = List.empty_list(unicode_type)
            for const in constituents:  function.append(tracer["parameters"]["metabolic_release"]["function"][const])

            self.metabolism_ids = List(["constituents","destination","function"])
            self.metabolism_params = List.empty_list(types.ListType(unicode_type))
            self.metabolism_params.append(constituents)
            self.metabolism_params.append(destination)
            self.metabolism_params.append(function)

        # Mortality
        if "mortality" in tracer["parameters"]:
            # Create list of mortality rates
            mort_rate = []  # [linear,quadratic,oxygen]
            if "linear" in tracer["parameters"]["mortality"]["mortality_rate"]:     mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["linear"]))
            else:   mort_rate.append(np.float64(0.))
            if "quadratic" in tracer["parameters"]["mortality"]["mortality_rate"]:  mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["quadratic"]))
            else:   mort_rate.append(np.float64(0.))
            if "oxygen" in tracer["parameters"]["mortality"]["mortality_rate"]:     mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["oxygen"]))
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
            
        # # Mortality
        # if "mortality" in tracer["parameters"]:
        #     # Create float option numbers for use in numba typed.List
        #     if "oxygen_limited" in tracer["parameters"]["mortality"]:
        #         # [0] False, [1] True
        #         if tracer["parameters"]["mortality"]["oxygen_limited"]:     tracer["parameters"]["mortality"]["oxygen_limited"] = 1
        #         else:   tracer["parameters"]["mortality"]["oxygen_limited"] = 0
        #     else:   tracer["parameters"]["mortality"]["oxygen_limited"] = 0     # Default to False if not in parameter list

        #     if "temp_limited" in tracer["parameters"]["mortality"]:
        #         # [0] False, [1] True
        #         if tracer["parameters"]["mortality"]["temp_limited"]:       tracer["parameters"]["mortality"]["temp_limited"] = 1
        #         else:   tracer["parameters"]["mortality"]["temp_limited"] = 0
        #     else:   tracer["parameters"]["mortality"]["temp_limited"] = 0     # Default to False if not in parameter list

        #     # Create list of mortality rates
        #     mort_rate = []  # [linear,quadratic]
        #     if "linear" in tracer["parameters"]["mortality"]["mortality_rate"]: mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["linear"]))
        #     else:   mort_rate.append(np.float64(0.))
        #     if "quadratic" in tracer["parameters"]["mortality"]["mortality_rate"]: mort_rate.append(np.float64(tracer["parameters"]["mortality"]["mortality_rate"]["quadratic"]))
        #     else:   mort_rate.append(np.float64(0.))

        #     tracer["parameters"]["mortality"]["mortality_rate"] = np.array(mort_rate,dtype=np.float64)

        #     self.mortality_ids = List.empty_list(unicode_type)
        #     self.mortality_params = List.empty_list(float64[:])
            
        #     for key,val in tracer["parameters"]["mortality"].items():
        #         self.mortality_ids.append(key)
        #         if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
        #         self.mortality_params.append(val)

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
            self.respiration_params = List.empty_list(float64)

            for key,val in tracer["parameters"]["respiration"].items():
                self.respiration_ids.append(key)
                if not isinstance(val, np.float64): val = np.float64(val)   # Convert type to array of floats for typed.List
                self.respiration_params.append(val)

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
        else:
            dict_type = types.DictType(types.unicode_type, types.float64)
            self.om_partition = Dict.empty(key_type=unicode_type,value_type=dict_type)

            # self.om_partition = List.empty_list(dict_type)
            # self.om_partition = lst

        # Oxygen inhibition
        if "oxygen_inhibition" in tracer["parameters"]:
            self.oxygen_limited = tracer["parameters"]["oxygen_inhibition"]["oxygen_limited"]
            # Initialize oxygen inhibition factor to array of 1. (if phytoplankton is temperature limited this will be updated later otherwise will stay as 1.)
            if num_layers > 1:  self.oxy_limitation_factor = np.ones(num_layers-1, dtype=np.float64)
            else:   self.oxy_limitation_factor = np.float64(1.)

            if "half_sat" in tracer["parameters"]["oxygen_inhibition"]: self.half_sat_oxygen = tracer["parameters"]["oxygen_inhibition"]["half_sat"]
            else:   self.half_sat_oxygen = 0.

            if "exponent" in tracer["parameters"]["oxygen_inhibition"]: self.oxy_limitation_exponent = tracer["parameters"]["oxygen_inhibition"]["exponent"]
            else:   self.oxy_limitation_exponent = 1.

            # if "oxygen_limited" in tracer["parameters"]["oxygen_inhibition"]:   self.oxygen_limited = tracer["parameters"]["oxygen_inhibition"]["oxygen_limited"]
            # else:   self.oxygen_limited = False

            # if self.oxygen_limited:
            #     # Default Hill exponent to 1 if not included in parameter list
            #     if "exponent" not in tracer["parameters"]["oxygen_inhibition"]:     tracer["parameters"]["oxygen_inhibition"]["exponent"] = np.float64(1.)

            #     self.oxy_inhib_ids = List.empty_list(unicode_type)
            #     self.oxy_inhib_params = List.empty_list(float64[:])
                
            #     for key,val in tracer["parameters"]["oxygen_inhibition"].items():
            #         self.oxy_inhib_ids.append(key)
            #         if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
            #         self.oxy_inhib_params.append(val)
        else:   
            self.oxygen_limited = False
            self.oxy_limitation_factor = 1.

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
        else:   
            self.temp_limited = False
            self.temp_regulation_factor = 1.

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
        
        # Reorder reactions (grazing needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "respiration"] + [item for item in self.reactions if item["type"] != "respiration"]
        self.reactions = [item for item in self.reactions if item["type"] == "grazing"] + [item for item in self.reactions if item["type"] != "grazing"]
        
        # Boolean to determine whether activity and basal respiration will be calculated for zooplankton processes
        self.calc_respiration = False
        for reac in self.reactions:
            if reac["type"] == "respiration":
                self.calc_respiration = True
                break
        # Boolean to determine whether grazing rates will be calculated for zooplankton processes
        self.calc_grazing = False
        for reac in self.reactions:
            if reac["type"] == "grazing":
                self.calc_grazing = True
                break
  

    def zoo(self, iter, base_element, temperature, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers):
        
        # Zero out grazing rates and total ingestion
        for prey in self.grazing_rates:
            self.grazing_rates[prey] = np.zeros_like(conc[tracer_map[self.abbrev][0]])
        self.total_ingestion = np.zeros_like(self.conc_ratio, dtype=np.float64)
        
        # Flag to return total ingestion
        return_ingestion = True     # True if no grazing rates have been calculated

        # Calculate temp regulation factor (if necessary)
        if self.temp_limited:
            # self.temp_regulation_factor = temperature_dependence(base_temp, temperature, self)
            self.temp_regulation_factor = temperature_dependence(temperature, self.temp_reg_ids, self.temp_reg_params)

        # Calculate oxygen limitation factor (if necessary)
        if self.oxygen_limited:
            self.oxy_limitation_factor = monod(conc[tracer_map["o2"][0]], self.half_sat_oxygen, self.oxy_limitation_exponent)
        
        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "grazing" and self.abbrev in reac["produced"].keys():       
                prey_composition = tracers[c[0]].composition
                zoo_composition = self.composition
                zoo_index = list(p).index(self.abbrev)
                if len(p) > 1:  # if prey constituent doesn't exist in zooplankton, attribute it to organic matter pool
                    if zoo_index == 0:  om_index = 1
                    else:               om_index = 0
                    om_composition = tracers[p[om_index]].composition
                else:   om_composition = List.empty_list(unicode_type)
                if return_ingestion:
                    self.total_ingestion = self.grazing(c, p, ec, ep, ic, ip, self.grazing_ids, self.grazing_params, self.grazing_preferences, self.grazing_rates, self.prey_availability, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, zoo_index, om_composition, prey_composition, zoo_composition, ingestion_flag = 1)
                    return_ingestion = False    # Turn off flag so total ingestion isn't calculated again (reduces overhead cost if there are multiple grazing reactions)
                else:   self.grazing(c, p, ec, ep, ic, ip, self.grazing_ids, self.grazing_params, self.grazing_preferences, self.grazing_rates, self.prey_availability, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, zoo_index, om_composition, prey_composition, zoo_composition, ingestion_flag = 0)
            if reac["type"] == "egestion":
                om_composition = tracers[p[0]].composition
                zoo_composition = self.composition
                # self.egestion(base_element, c, p, ec, ep, self.assimilation_efficiency, self.ingestion_efficiency, self.total_ingestion, self.om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition)
                self.egestion(base_element, c, p, ec, ep, self.unassimilated_egestion, self.assimilation_efficiency, self.ingestion_efficiency, self.total_ingestion, self.om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition)
            if reac["type"] == "excretion":     
                # if self.calc_respiration:   self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basal_respiration)
                # else:                       self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, 0., 0.)
                if self.calc_respiration:   pass
                else:
                    activity_respiration = 0.
                    basal_respiration = 0.
                self.excretion(base_element, c, p, ec, ep, self.cell_quota_ids, self.cell_quota_opt, self.excretion_ids, self.excretion_params, activity_respiration, basal_respiration, self.ingestion_efficiency, self.total_ingestion, conc, conc_ratio, d_dt, tracer_map, self.composition)
            if reac["type"] == "metabolic_release":
                self.metabolic_release(c, p, ec, ep, self.cell_quota_ids, self.cell_quota_opt, self.metabolism_ids, self.metabolism_params, self.basal_metabolic_rate, self.assimilation_efficiency, self.ingestion_efficiency, self.total_ingestion, self.oxy_limitation_factor, self.temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, self.composition)
            if reac["type"] == "mortality":     
                om_composition = tracers[p[0]].composition
                zoo_composition = self.composition
                self.mortality(c, p, ec, ep, ic, ip, self.mortality_ids, self.mortality_params, self.oxy_limitation_factor, self.temp_regulation_factor, self.om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition)
            if reac["type"] == "respiration":   activity_respiration, basal_respiration = self.respiration(self.abbrev, base_element, c, p, self.respiration_ids, self.respiration_params, self.assimilation_efficiency, self.ingestion_efficiency, self.total_ingestion, self.temp_regulation_factor, conc, d_dt, tracer_map, self.composition)
        

    def add_prey(self, prey):
        for p in prey:
            self.grazing_rates[p] = np.zeros(self.conc.shape[1], dtype=np.float64)
            self.prey_availability[p] = np.zeros(self.conc.shape[1], dtype=np.float64)


    @staticmethod
    @njit
    def egestion(base_element, c, p, ec, ep, unassimilated_egestion, assimilation_efficiency, ingestion_efficiency, total_ingestion, om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition):
    # def egestion(base_element, c, p, ec, ep, assimilation_efficiency, ingestion_efficiency, total_ingestion, om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition):
        """
        Definition:: Calculates zooplankton loss to detrital pool as feacal pellet production
        Return:: Egestion rate
        """
        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]

        # Get index of base element in zooplankton
        # base_index = zoo_composition.index(base_element)

        # Calculate egestion rate
        egestion = ingestion_efficiency * total_ingestion
        if len(unassimilated_egestion) > 0:
            for const in unassimilated_egestion:
                index = zoo_composition.index(const)
                egestion[index] *= (1. - assimilation_efficiency)
        # if egestion_formulation == "split": # egestion not fully assimilated
        #     egestion[base_index] = egestion[base_index] * ( 1. - assimilation_efficiency )

        # Calculate concentration ratios
        ratios = np.zeros((len(zoo_composition),len(egestion[0])),dtype=np.float64)
        for const in zoo_composition:
            if const in om_composition:
                index_zoo = zoo_composition.index(const)
                index_om = om_composition.index(const)
                ratios[index_om] = conc_ratio[tracer_map[cons][index_zoo]]

        # Update d_dt
        if prod in om_partition:    # Apply partition (if necessary)
            for const in om_partition[prod]:
                om_const_index = om_composition.index(const)    # Get index of constituent in organic matter pool
                if const in zoo_composition:    # Only apply rate if constituent is also in zooplankton
                    zoo_const_index = zoo_composition.index(const)  # Get index of constituent in zooplankton
                    # d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * conc_ratio[tracer_map[cons][zoo_const_index]] * egestion[zoo_const_index] * om_partition[prod][const]
                    # # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * conc_ratio[tracer_map[prod][om_const_index]] * egestion[zoo_const_index] * om_partition[prod][const]
                    # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * egestion[zoo_const_index] * om_partition[prod][const]
                    d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * egestion[zoo_const_index] * om_partition[prod][const]
                    # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * conc_ratio[tracer_map[prod][om_const_index]] * egestion[zoo_const_index] * om_partition[prod][const]
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * egestion[zoo_const_index] * om_partition[prod][const]

        else:
            for const in om_composition:
                om_const_index = om_composition.index(const)
                if const in zoo_composition:
                    zoo_const_index = zoo_composition.index(const)
                    # d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * conc_ratio[tracer_map[cons][zoo_const_index]] * egestion[zoo_const_index]
                    # # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * conc_ratio[tracer_map[prod][om_const_index]] * egestion[zoo_const_index]
                    # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * egestion[zoo_const_index]
                    d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * egestion[zoo_const_index]
                    # d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * conc_ratio[tracer_map[prod][om_const_index]] * egestion[zoo_const_index]
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * egestion[zoo_const_index]

    
    @staticmethod
    @njit
    def excretion(base_element, c, p, ec, ep, cell_quota_ids, cell_quota_opt, excretion_ids, excretion_params, activity_respiration, basal_respiration, ingestion_efficiency, total_ingestion, conc, conc_ratio, d_dt, tracer_map, composition):
        """
        Definition:: Calculates excretion of zooplankton to nutrient or detrital pool.
                     Excretion can be represented either as a constant rate or as a fraction of zooplankton grazing on phytoplankton.
        Return:: Excretion rate
        """
        # Extract parameters
        constituents = excretion_params[excretion_ids.index("constituents")]
        destination = excretion_params[excretion_ids.index("destination")]
        excretion_function = excretion_params[excretion_ids.index("function")]

        # # Nutrient elements
        # nutrient_elements = Dict.empty(key_type=types.unicode_type,value_type=types.unicode_type)
        # nutrient_elements["no3"] = "n"
        # nutrient_elements["nh4"] = "n"
        # nutrient_elements["po4"] = "p"

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_p = ep[prod]
        ind_p = list(elem_p).index(1.)

        # Get the element of the destination tracer
        element_index = destination.index(prod)
        element = constituents[element_index]
        
        # Get concentration of constituent in zooplankton
        # nutrient_index = composition.index(nutrient_elements[prod])     # Index of excreted nutrient
        # zoo = conc[tracer_map[cons][0]]
        nutrient_index = composition.index(element)         # Index of excreted nutrient
        zoo = conc[tracer_map[cons][nutrient_index]]
        base_index = composition.index(base_element)    # Index of base element

        # Extract cell quota id
        # quota_index = cell_quota_ids.index(nutrient_elements[prod])
        quota_index = cell_quota_ids.index(element)

        # Extract parameter indices
        # nutrient_excretion = excretion_ids.index(prod)

        # Excretion rate based on grazing rate
        # if excretion_params[nutrient_excretion] == -1.:
        if excretion_function[element_index] == "activity":
                excreted_base = np.maximum(np.zeros_like(total_ingestion[nutrient_index]), total_ingestion[base_index] * (1. - ingestion_efficiency) - activity_respiration)
                excreted_nutrient = np.maximum(np.zeros_like(total_ingestion[nutrient_index]), ( total_ingestion[nutrient_index] * (1. - ingestion_efficiency) ) + ( basal_respiration * conc_ratio[tracer_map[cons][nutrient_index]] ))
            
                excretion = np.maximum(np.zeros_like(excreted_base), excreted_nutrient/(excreted_base + 1.E-20) - cell_quota_opt[quota_index] ) * excreted_base

        else:   # constant excretion rate
            rate_str = excretion_function[element_index]
            excretion_rate = string_to_float(rate_str)
            # excretion_rate = np.float64(function[element_index])
            excretion = excretion_rate * zoo

            # Calculate excretion of excess nutrient if necessary
            if element != base_element:
                element_ratio = conc_ratio[tracer_map[cons][nutrient_index]]
                excretion = excretion * np.maximum(0., element_ratio - cell_quota_opt[quota_index])

            # excretion = excretion_params[nutrient_excretion] * zoo

            # # Calculate excretion of excess nutrient
            # element_ratio = zoo / conc[tracer_map[cons][base_index]]
            # excretion = excretion * np.maximum(0., element_ratio - cell_quota_opt[quota_index])

        # Update d_dt
        d_dt[tracer_map[cons][nutrient_index]] -= excretion
        d_dt[tracer_map[prod][ind_p]] += excretion


    @staticmethod
    @njit
    def grazing(c, p, ec, ep, ic, ip, grazing_ids, grazing_params, grazing_preferences, grazing_rates, prey_availability, temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, zoo_index, om_composition, prey_composition, zoo_composition, ingestion_flag):
        """
        Definition:: Calculates the zooplankton grazing rate on a particular species using user choice of the Ivlev Equation,
                    Holling Type I Response, Holling Type II Response, or Holling Type III Response.
        Return:: grazing_c - Grazing rate to be added to "grazing_rates[]" for later summation
                grazing_p - Porition of grazing rate allocated to zooplaknton (scaled by assimilation and ingestion efficiencies).
        """
        # Extract parameter indices
        max_grazing_rate = grazing_ids.index("max_grazing_rate")
        feeding_model = grazing_ids.index("feeding_model")
        if grazing_params[feeding_model] == 1.:     search_volume = grazing_ids.index("search_volume")
        else:   half_sat_grazing = grazing_ids.index("half_sat_grazing") 

        function = grazing_ids.index("function")    # grazing function
        if function != 4.:  # function not ivlev --> holling type grazing functions
            if grazing_params[feeding_model] == 1:  # clearance_rate based feeding model does not use feeding threshold
                feeding_threshold = 1.E-20  # Small constant to prevent divide by 0.
            else:
                feeding_threshold_id = grazing_ids.index("feeding_threshold")
                feeding_threshold = grazing_params[feeding_threshold_id]   # extract parameter to use for capture efficiency
        else:   # ivlev
            ivlev = grazing_ids.index("ivlev")
            feeding_threshold = 1.E-20  # Small constant to prevent divide by 0.

        # Extract dict
        cons = c[0]
        if len(p) > 1:  # if prey constituent doesn't exist in zooplankton, attribute it to organic matter pool
            prod = p[zoo_index]
        else:
            prod = p[0]

        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons][0]
        ind_p = ip[prod][0]

        tc = conc[tracer_map[cons][ind_c]]
        tp = conc[tracer_map[prod][ind_p]]

        # Calculate total food availability
        for prey, preference in grazing_preferences.items():
            # Concentration of base element in prey
            conc_prey = conc[tracer_map[prey][ind_c]]

            # Capture efficiency for current prey in list of available
            # eff_prey = conc_prey / ( conc_prey + grazing_params[feeding_threshold] )
            if grazing_params[feeding_model] == 2.:     # "half_saturation" feeding method uses capture efficiency to scale prey availability
                eff_prey = conc_prey / ( conc_prey + feeding_threshold )    # Changed to this after realizing code would break if using ivlev grazing function
            # else:   eff_prey = 1.   # "clearance_rate" feeding method does not use captre efficiency (set to one for no scaling)
            else:   eff_prey = np.ones_like(conc_prey)   # "clearance_rate" feeding method does not use captre efficiency (set to one for no scaling)

            # Total food availability is sum for all prey (prey concentration squared for ONLY holling-3 sigmoidal behavior)
            if grazing_params[function] == 3.:  prey_availability[prey] = preference * eff_prey * (conc_prey**2)
            else:   prey_availability[prey] = preference * eff_prey * conc_prey

        # Total food availability is sum for all prey
        total_available = np.zeros_like(tc)
        for prey in prey_availability:  total_available += prey_availability[prey]

        # Calculate grazing function
        if grazing_params[function] == 1.:    # holling-1, Linear
            # Calculate slope based on prey concentration 
            slope = np.zeros_like(tc)
            if grazing_params[feeding_model] == 1:  # clearance_rate
                for i in range(0,len(slope)):
                    if tc[i] < (2 / grazing_params[search_volume]):    slope[i] = (grazing_params[max_grazing_rate] * grazing_params[search_volume]) / 2
                    else:   slope[i] = grazing_params[max_grazing_rate]
            else:   # half_saturation
                for i in range(0,len(slope)):
                    if tc[i] < (2 * grazing_params[half_sat_grazing]):    slope[i] = grazing_params[max_grazing_rate]/(2 * grazing_params[half_sat_grazing])
                    else:   slope[i] = grazing_params[max_grazing_rate]
            
            # Calculate specific grazing rate for individual prey
            grazing = slope * prey_availability[cons] * tp

            # Calculate total uptake rate
            total_uptake = slope * total_available * tp

        elif grazing_params[function] == 2.:    # holling-2, Hyperbolic
            # Calculate specific grazing rate for individual prey
            # grazing = ( grazing_params[max_grazing_rate] * prey_availability[cons] ) / ( total_available + grazing_params[half_sat_grazing] ) * tp

            # Calculate total uptake rate
            if grazing_params[feeding_model] == 1:  # clearance_rate
                # grazing = ( grazing_params[max_grazing_rate] * grazing_params[search_volume] * prey_availability[cons] ) / ( (grazing_params[search_volume] * total_available) + grazing_params[max_grazing_rate] ) * tp
                # total_uptake = ( grazing_params[max_grazing_rate] * grazing_params[search_volume] * total_available ) / ( (grazing_params[search_volume] * total_available) + grazing_params[max_grazing_rate] ) * tp

                total_uptake = ( grazing_params[max_grazing_rate] * grazing_params[search_volume] * total_available ) / ( (grazing_params[search_volume] * total_available) + grazing_params[max_grazing_rate] ) * tp
                grazing = ( total_uptake * prey_availability[cons] ) / ( 1.E-20 + total_available )
            else:   # half_saturation
                grazing = ( grazing_params[max_grazing_rate] * prey_availability[cons] ) / ( total_available + grazing_params[half_sat_grazing] ) * tp
                total_uptake = ( grazing_params[max_grazing_rate] * total_available ) / ( total_available + grazing_params[half_sat_grazing] ) * tp
                        
        elif grazing_params[function] == 3.:    # holling-3, Sigmoidal
            # # Calculate specific grazing rate for individual prey, half saturation constant is also squared for sigmoidal behavior (squared prey concentration handled in calculation of prey availability)
            # grazing = ( grazing_params[max_grazing_rate] * prey_availability[cons] ) / ( total_available + (grazing_params[half_sat_grazing]**2) ) * tp

            # # Calculate total uptake rate
            # total_uptake = ( grazing_params[max_grazing_rate] * total_available ) / ( (total_available**2) + (grazing_params[half_sat_grazing]**2) ) * tp

            if grazing_params[feeding_model] == 1:  # clearance_rate
                grazing = ( grazing_params[max_grazing_rate] * grazing_params[search_volume] * prey_availability[cons] ) / ( (grazing_params[search_volume] * total_available) + (grazing_params[max_grazing_rate]**2) ) * tp
                total_uptake = ( grazing_params[max_grazing_rate] * grazing_params[search_volume] * total_available ) / ( (grazing_params[search_volume] * total_available)**2 + (grazing_params[max_grazing_rate]**2) ) * tp
            else:   # half_saturation
                grazing = ( grazing_params[max_grazing_rate] * prey_availability[cons] ) / ( total_available + (grazing_params[half_sat_grazing]**2) ) * tp
                total_uptake = ( grazing_params[max_grazing_rate] * total_available ) / ( (total_available**2) + (grazing_params[half_sat_grazing]**2) ) * tp

        elif grazing_params[function] == 4.:      # ivlev, Exponential
            # Calculate specific grazing rate for individual prey
            grazing = grazing_params[max_grazing_rate] * ( 1 - np.exp( -grazing_params[ivlev] * total_available ) ) * ( prey_availability[cons]/total_available ) * tp

            # Calculate total uptake rate (total_available/total_available = 1 so unnecessary multiplication left off)
            total_uptake = grazing_params[max_grazing_rate] * ( 1 - np.exp( -grazing_params[ivlev] * total_available ) )  * tp
                
        # Temperature regulation
        grazing = grazing * temp_regulation_factor
        total_uptake = total_uptake * temp_regulation_factor

        # Extract cell quotas from prey
        ratios = np.zeros((len(zoo_composition),len(tp)),dtype=np.float64)
        for const in zoo_composition:
            if const in prey_composition:
                index_prey = prey_composition.index(const)
                index_pred = zoo_composition.index(const)
                ratios[index_pred] = conc_ratio[tracer_map[cons][index_prey]]

        # Update d_dt
        for i in range(len(elem_c)):    d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * grazing

        for j in range(len(elem_p)):    
            d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * grazing

            # Update grazing rate dictionary
            grazing_rates[cons] = elem_p[j] * ratios[j] * grazing
        
        # if prey constituent doesn't exist in zooplankton, attribute it to organic matter pool
        if len(p) > 1:
            om = p[1-zoo_index] 
            for const in prey_composition:
                if const not in zoo_composition and const in om_composition:    # Check if constituent in om_composition and not in zooplankton
                    const_index_prey = prey_composition.index(const)    # constituent index in prey
                    const_index_om = om_composition.index(const)        # constituent index in organic matter

                    # Apply rate to organic matter
                    d_dt[tracer_map[om][const_index_om]] += conc_ratio[tracer_map[cons][const_index_prey]] * grazing

        # Calculate total ingestion (if necessary)
        if ingestion_flag == 1:     # This is the first time the function has been called during this time step, calculate and return total ingestion
            # Create array for total ingestion of all prey 
            # NEED to calculate total instead of specific because cannibalism may not be included in reaction list (will lead to error in future reactions if cannibalism is excluded)
            total_ingestion = np.zeros((len(zoo_composition),len(tp)),dtype=np.float64)
            for prey in prey_availability:
                for const in prey_composition:
                    if const in zoo_composition:
                        index_prey = prey_composition.index(const)
                        index_pred = zoo_composition.index(const)
                        # total_ingestion[index_pred] += temp_regulation_factor * (total_uptake / (total_available + 1.E-20)) * conc_ratio[tracer_map[prey][index_prey]] * prey_availability[prey]
                        total_ingestion[index_pred] += (total_uptake / (total_available + 1.E-20)) * conc_ratio[tracer_map[prey][index_prey]] * prey_availability[prey]
        
            return total_ingestion


    @staticmethod
    @njit
    def metabolic_release(c, p, ec, ep, cell_quota_ids, cell_quota_opt, metabolism_ids, metabolism_params, basal_metabolic_rate, assimilation_efficiency, ingestion_efficiency, total_ingestion, oxy_limitation_factor, temp_regulation_factor, conc, conc_ratio, d_dt, tracer_map, composition):
        """
        Definition:: Calculates metabolic release of excess non-limiting nutrients of zooplankton to nutrient or detrital pool.
        Return:: 
        """
        # Extract parameters
        constituents = metabolism_params[metabolism_ids.index("constituents")]
        destination = metabolism_params[metabolism_ids.index("destination")]
        function = metabolism_params[metabolism_ids.index("function")]

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = list(elem_c).index(1.)
        ind_p = list(elem_p).index(1.)

        # Get the element of the destination tracer
        element_index = destination.index(prod)
        element = constituents[element_index]

        # Get concentration of element in zooplankton
        nutrient_index = composition.index(element)         # Index of excreted nutrient
        zoo = conc[tracer_map[cons][nutrient_index]]

        # Calculate energy cost of ingestion
        energy_cost = 1. - assimilation_efficiency - ingestion_efficiency

        # Check assimilation rate
        assimilation_rate = np.zeros((len(composition),len(zoo)),dtype=np.float64)
        for i in range(len(composition)):
            if i == 0:  # index of base element
                assimilation_rate[i] = assimilation_efficiency * total_ingestion[i]
            else:   # include energy cost of ingestion for non-base element constituents
                assimilation_rate[i] = ( assimilation_efficiency + energy_cost ) * total_ingestion[i]

        # Calculate concentration ratios in assimilation rate
        assimilation_ratios = np.zeros((len(composition),len(zoo)),dtype=np.float64)
        for i in range(len(composition)):
            assimilation_ratios[i] = assimilation_rate[i] / ( conc_ratio[tracer_map[cons][i]] + 1.E-20 )

        # Set assimilation ratio of base element index to ~0 for limiting element determination
        assimilation_ratios[0,:] = 1.E-20

        # Find index of limiting nutrient at all depths
        limiting_index = np.argmin(assimilation_ratios,axis=0)
        correction = np.zeros((len(composition),len(zoo)),dtype=np.float64)

        for depth in range(0,len(zoo)):
            limiting_element_index = limiting_index[depth]

            if limiting_element_index == 0:  # base element is limiting element
                for i in range(len(composition)):
                    if i > 0:   # only perform calculation for non-base element indices
                        quota_index = cell_quota_ids.index(composition[i])
                        correction[i,depth] = max(0., (1. - ingestion_efficiency) * total_ingestion[i,depth] - (cell_quota_opt[quota_index] * assimilation_rate[limiting_element_index,depth]) )

            elif ( limiting_element_index != 0 ) and ( assimilation_rate[limiting_element_index,depth] < conc_ratio[tracer_map[cons][limiting_element_index],depth] ):
                for i in range(len(composition)):
                    if i == 0:  # perform caculation for base element first
                        quota_index = cell_quota_ids.index(composition[limiting_element_index])
                        correction[i,depth] = max(0., assimilation_rate[i,depth] - ( (1. - ingestion_efficiency) * total_ingestion[limiting_element_index,depth] / cell_quota_opt[quota_index] ) )
                    elif ( i > 0 ) and ( i != limiting_element_index ):  # perform calculation for other non-limiting elements
                        quota_index = cell_quota_ids.index(composition[i])
                        correction[i,depth] = max(0., (1. - ingestion_efficiency) * total_ingestion[i,depth] - cell_quota_opt[quota_index] * (assimilation_rate[0,depth] - correction[0,depth]))
        
        # Calculate metabolic release
        if function[element_index] == "activity":
            release = basal_metabolic_rate * oxy_limitation_factor * temp_regulation_factor * zoo + correction[nutrient_index]

        elif function[element_index] == "excess_only":
            release = correction[nutrient_index]

        # Update d_dt
        d_dt[tracer_map[cons][nutrient_index]] -= release
        d_dt[tracer_map[prod][ind_p]] += release


    @staticmethod
    @njit
    def mortality(c, p, ec, ep, ic, ip, mortality_ids, mortality_params, oxy_limitation_factor, temp_regulation_factor, om_partition, conc, conc_ratio, d_dt, tracer_map, om_composition, zoo_composition):
        """
        Definition:: Calculates the non-grazing mortality of planktoninc species
        """
        # Extract parameter indices
        mortality_rate = mortality_ids.index("mortality_rate")
        # oxygen_limited = mortality_ids.index("oxygen_limited")
        temp_limitation = mortality_ids.index("temp_limitation")
        # if mortality_params[oxygen_limited][0] == 1.:   # include an additional oxygen-dependent mortality term
        #     mortality_rate_oxy = mortality_ids.index("mortality_rate_oxy")
        #     # half_sat_oxygen = mortality_ids.index("half_sat_oxygen")

        # Extract dict
        cons = c[0]
        prod = p[0]
        elem_c = ec[cons]
        elem_p = ep[prod]
        ind_c = ic[cons][0]
        ind_p = ip[prod][0]
        
        zoo = conc[tracer_map[cons][ind_c]]

        # Calculate mortality rate
        # mortality = ( mortality_params[mortality_rate][0] * zoo ) + ( mortality_params[mortality_rate][1] * (zoo**2) )

        linear = mortality_params[mortality_rate][0] * zoo
        quadratic = mortality_params[mortality_rate][1] * (zoo**2)
        oxygen = ( 1. - oxy_limitation_factor ) * mortality_params[mortality_rate][2] * zoo

        # Apply temperature limitation if necessary (1. = True)
        if mortality_params[temp_limitation][0] == 1.:  linear *= temp_regulation_factor
        if mortality_params[temp_limitation][1] == 1.:  quadratic *= temp_regulation_factor
        if mortality_params[temp_limitation][2] == 1.:  oxygen *= temp_regulation_factor

        mortality = linear + quadratic + oxygen
        
        # # Oxygen limitation
        # if mortality_params[oxygen_limited][0] == 1.:
        #     # oxy_limitation_factor = np.minimum(1., nutrient_limitation(conc[tracer_map["o2"][0]], mortality_params[half_sat_oxygen]))
        #     mortality += (1. - oxy_limitation_factor) * mortality_params[mortality_rate_oxy] * zoo
    
        # Calculate concentration ratios
        ratios = np.zeros((len(zoo_composition),len(zoo)),dtype=np.float64)
        for const in zoo_composition:
            if const in om_composition:
                index_zoo = zoo_composition.index(const)
                index_om = om_composition.index(const)
                ratios[index_om] = conc_ratio[tracer_map[cons][index_zoo]]

        # Update d_dt
        if prod in om_partition:    # Apply partition (if necessary)
            for const in om_partition[prod]:
                om_const_index = om_composition.index(const)    # Get index of constituent in organic matter pool
                if const in zoo_composition:    # Only apply rate if constituent is also in zooplankton
                    zoo_const_index = zoo_composition.index(const)  # Get index of constituent in zooplankton
                    d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * conc_ratio[tracer_map[cons][zoo_const_index]] * mortality * om_partition[prod][const]
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * mortality * om_partition[prod][const]

        else:
            for const in om_composition:
                om_const_index = om_composition.index(const)
                if const in zoo_composition:
                    zoo_const_index = zoo_composition.index(const)
                    d_dt[tracer_map[cons][zoo_const_index]] -= elem_c[zoo_const_index] * conc_ratio[tracer_map[cons][zoo_const_index]] * mortality
                    d_dt[tracer_map[prod][om_const_index]] += elem_p[om_const_index] * ratios[om_const_index] * mortality


    @staticmethod
    @njit
    def respiration(abbrev, base_element, c, p, respiration_ids, respiration_params, assimilation_efficiency, ingestion_efficiency, total_ingestion, temp_regulation_factor, conc, d_dt, tracer_map, composition):
        """
        Definition:: Calculates zooplankton respiration
        """

        # Extract parameter indices
        respiration_rate = respiration_ids.index("respiration_rate")
        conv_o2 = False
        if "convert_o2" in respiration_ids:
            conv_o2 = True
            convert_o2 = respiration_ids.index("convert_o2")
        conv_co2 = False
        if "convert_co2" in respiration_ids:
            conv_co2 = True
            convert_co2 = respiration_ids.index("convert_co2")

        # Locate index of base element
        base_index = composition.index(base_element)

        # Get concentration of base element
        zoo = conc[tracer_map[abbrev][base_index]]
        
        # Calculate respiration rates
        activity_respiration = (1 - assimilation_efficiency - ingestion_efficiency) * total_ingestion[base_index]
        basal_respiration = temp_regulation_factor * respiration_params[respiration_rate] * zoo

        # Update d_dt
        total_respiration = activity_respiration + basal_respiration
        d_dt[tracer_map[abbrev][base_index]] -= total_respiration
        if c and "o2" in c:   
            if conv_o2:     d_dt[tracer_map["o2"][0]] -= total_respiration * respiration_params[convert_o2]
            else:           d_dt[tracer_map["o2"][0]] -= total_respiration

        if p and "co2" in p:
            if conv_co2:    d_dt[tracer_map["co2"][0]] += total_respiration * respiration_params[convert_co2]
            else:           d_dt[tracer_map["co2"][0]] += total_respiration

        return activity_respiration, basal_respiration
    