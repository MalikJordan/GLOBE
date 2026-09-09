import os
import sys
import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from fractions import Fraction
from functions.seasonal_cycling import *
from functions.other_functions import tracer_elements, temperature_dependence, monod
from fractions import Fraction
np.set_printoptions(precision=20)
class Detritus():
    """
    
    """

    def __init__(self, abbrev, base_element, physical, reactions, **tracer):
        
        # Variales that will be used later ---------------------------------------------------------------
        num_layers = physical["water_column"]["num_layers"]
        iters = physical["simulation"]["iters"]
        composition = physical["initial_conditions"][abbrev]["composition"]     # Initial concentrations
        if "scale" in physical["initial_conditions"][abbrev]:   scale = physical["initial_conditions"][abbrev]["scale"]     # Scaling factor (if initial concentration is split between multiple detrital groups)
        else:   scale = 1.  # No scaling

        # Add important keys ---------------------------------------------------------------
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]
        self.form = tracer["form"]

        # Add parameters ---------------------------------------------------------------
        # Light attenuation
        if "light_attenuation" in tracer["parameters"]:
            self.light_attenuation = np.array([tracer["parameters"]["light_attenuation"]],dtype=np.float64)
        else:
            self.light_attenuation = np.array([0.],dtype=np.float64)

        # Dissolution
        if "dissolution" in tracer["parameters"]:
            self.dissolution_ids = List.empty_list(unicode_type)
            self.dissolution_params = []

            for key,val in tracer["parameters"]["dissolution"]["dissolution_rate"].items():
                self.dissolution_ids.append(key)
                self.dissolution_params.append(np.float64(tracer["parameters"]["dissolution"]["dissolution_rate"][key]))

        # Loss (Generic)
        if "loss" in tracer["parameters"]:
            loss_dict_type = types.DictType(types.unicode_type, types.float64)
            self.loss_ids = List.empty_list(unicode_type)   # stores produced tracer names (loss parameters will be saved for each tracer individually)
            self.loss_params = List.empty_list(loss_dict_type)

            for outer_key,inner_dict in tracer["parameters"]["loss"].items():
                # Create temporary dictionary
                temp = Dict.empty(key_type=unicode_type, value_type=float64)

                # Add "exponent" to inner_dict (if necessary)
                if inner_dict["function"] == "half_saturation" and "exponent" not in inner_dict:    inner_dict["exponent"] = 1.

                # Add inner values to temporary dictionary
                for inner_key,inner_val in inner_dict.items():
                    # Create numeric codes for loss option string
                    if inner_key == "function":
                        if inner_val == "constant":             temp["function"] = 1.
                        elif inner_val == "half_saturation":    temp["function"] = 2.

                    # Other inner values are already ints/floats
                    else:
                        temp[inner_key] = np.float64(inner_val) # conversion to make sure all values are floats

                # Add outer_key to loss_ids and temp to loss_params
                self.loss_ids.append(outer_key)
                self.loss_params.append(temp)


        # Remineralization
        if "remineralization" in tracer["parameters"]:
            self.remineralization_ids = List.empty_list(unicode_type)
            self.remineralization_params = []

            if "convert_o2" in tracer["parameters"]["remineralization"]:
                # If conversion given as a Fraction string (ex: '1/2' instead of 0.5), convert to float
                if isinstance(tracer["parameters"]["remineralization"]["convert_o2"],str):
                    tracer["parameters"]["remineralization"]["convert_o2"] = np.float64(Fraction(tracer["parameters"]["remineralization"]["convert_o2"]))

            for key,val in tracer["parameters"]["remineralization"].items():
                self.remineralization_ids.append(key)
                self.remineralization_params.append(val)

        # Sedimentation
        if "sedimentation" in tracer["parameters"] and tracer["parameters"]["sedimentation"]["sinking"] == True:
            # Create attribute for background sinking velocity
            if num_layers > 1:
                self.sinking_velocity = np.ones(num_layers-1,dtype=np.float64) * tracer["parameters"]["sedimentation"]["background_sinking_rate"]
                self.sinking_velocity[-1] = np.float64(tracer["parameters"]["sedimentation"]["burial_velocity"])
            else:
                self.sinking_velocity = np.array([tracer["parameters"]["sedimentation"]["background_sinking_rate"]],dtype=np.float64)

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
        else:   self.temp_limited = False

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
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")

        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.initial_conc = np.zeros((len(self.composition),num_layers-1),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.initial_conc[const,:] = scale * conc[const][:-1] # Apply scaling factor here to prevent from applying multiple times in the above step
        else:   # Model as single box
            self.initial_conc = np.zeros((len(self.composition),num_layers),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.initial_conc[const,:] = scale * conc[const]      # Apply scaling factor here to prevent from applying multiple times in the above step
        
        # Add reactions ---------------------------------------------------------------
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

            # Delete "remineralization" reactions if this tracer is produced
            if ( reac["type"] == "remineralization" ) and ( abbrev in produced.keys() ):    self.reactions.pop()

            # Delete "grazing" reactions if this tracer is consumed
            if ( reac["type"] == "grazing" ) and ( abbrev in consumed.keys() ): self.reactions.pop()

            # Delete "uptake" reactions if this tracer is consumed
            if ( reac["type"] == "uptake" ) and ( abbrev in consumed.keys() ): self.reactions.pop()
        

    def detritus(self, base_element, temperature, conc, d_dt, tracer_map, tracers):

        # Calculate temp regulation factor (if necessary)
        if self.temp_limited:
            self.temp_regulation_factor = temperature_dependence(temperature, self.temp_reg_ids, self.temp_reg_params)
            
        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            
            if reac["type"] == "dissolution":       self.dissolution(c, p, ec, ep, ic, ip, self.dissolution_ids, self.dissolution_params, self.temp_regulation_factor, conc, d_dt, tracer_map, self.composition)
            if reac["type"] == "remineralization":  self.remineralization(c, p, ec, ep, ic, ip, self.remineralization_ids, self.remineralization_params, conc, d_dt, tracer_map)
    
    
    @staticmethod
    @njit
    def dissolution(c, p, ec, ep, ic, ip, dissolution_ids, dissolution_params, temp_regulation_factor, conc, d_dt, tracer_map, composition):

        # Extract dict
        cons = c[0]
        elem_c = ec[cons]
        ind_c = list(elem_c).index(1.)

        # Get concentration of nutreint in organic matter pool
        tc = conc[tracer_map[cons][ind_c]]

        prod = p[0]
        elem_p = ep[prod]

        # Extract parameter indices
        dissolution_rate_index = dissolution_ids.index(composition[ind_c])
        
        # Calculate dissolution rate
        dissolution = dissolution_params[dissolution_rate_index] * temp_regulation_factor * tc

        # Update d_dt
        for i in range(len(elem_c)):
            d_dt[tracer_map[cons][i]] -= elem_c[i] * dissolution
        for j in range(len(elem_p)):
            d_dt[tracer_map[prod][j]] += elem_p[j] * dissolution


    @staticmethod
    @njit
    def loss(c, p, ec, ep, ic, ip, loss_ids, loss_params, conc, conc_ratio, d_dt, tracer_map, cons_composition, prod_composition):
    
        # Source tracer of loss rate
        cons = c[0]         # Tracer
        elem_c = ec[cons]   # Affected constituents
        ind_c = ic[cons][0] # Base element index
    
        # Destination tracer of loss rate
        prod = p[0]         # Tracer
        elem_p = ep[prod]   # Affected constituents
        ind_p = ip[prod][0] # Base element index
    
        # Identifiy loss parameters for consumed tracer
        ids = loss_ids.index(prod)
        params = loss_params[ids]
    
        if params["function"] == 1.: # constant
            # Calculate loss rate
            loss = params["loss_rate"] * conc[tracer_map[cons][ind_c]]
    
        elif params["function"] == 2.: # half saturation
            # Calculate loss rate
            loss = params["loss_rate"] * monod(conc[tracer_map[cons][ind_c]], params["half_sat_loss"], params["exponent"]) * conc[tracer_map[cons][ind_c]]
    
        # Extract concentration ratios
        ratios = np.zeros((len(prod_composition),len(conc[tracer_map[prod][0]])),dtype=np.float64)
        for const in cons_composition:
            if const in prod_composition:
                index_cons = cons_composition.index(const)
                index_prod = prod_composition.index(const)
                ratios[index_prod] = conc_ratio[tracer_map[cons][index_cons]]
    
        # Update d_dt
        for i in range(0,len(elem_c)):
            d_dt[tracer_map[cons][i]] -= elem_c[i] * conc_ratio[tracer_map[cons][i]] * loss
        for j in range(0,len(elem_p)):
            d_dt[tracer_map[prod][j]] += elem_p[j] * ratios[j] * loss


    @staticmethod
    @njit
    def remineralization(c, p, ec, ep, ic, ip, remineralization_ids, remineralization_params, conc, d_dt, tracer_map):

        # Extract parameter indices
        remin_rate = remineralization_ids.index("remineralization_rate")
        convert = False
        if "convert_o2" in remineralization_ids:
            convert = True
            convert_o2 = remineralization_ids.index("convert_o2")
        
        # Extract dict
        if len(c) > 1 and "o2" in c:
        # Remineralization of carbon also affects oxygen (sink)
            for t in c:
                if t == "o2":   pass
                else:
                    cons = t
                    break
        else:   cons = c[0]
        elem_c = ec[cons]   # element consumed
        ind_c = list(elem_c).index(1.)  # index of element

        # Get concentration of remineralized nutrient in organic matter pool
        tc = conc[tracer_map[cons][ind_c]]

        if not p:    pass
        else:
            prod = p[0]
            elem_p = ep[prod]   # element produced
        
        # Calculate remineralization rate
        remineralization = (remineralization_params[remin_rate]) * tc

        # Update d_dt
        for i in range(len(elem_c)):
            d_dt[tracer_map[cons][i]] -= elem_c[i] * remineralization

        if "o2" in c:
            if convert: d_dt[tracer_map["o2"][0]] -= remineralization * remineralization_params[convert_o2]
            else:   d_dt[tracer_map["o2"][0]] -= remineralization
        
        if not p:    pass
        else:
            for j in range(len(elem_p)):
                d_dt[tracer_map[prod][j]] += elem_p[j] * remineralization
