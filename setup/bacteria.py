import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, monod, nutrient_limitation, tracer_elements
from fractions import Fraction
class Bacteria():
    """
    
    """


    def __init__(self, abbrev, base_element, iters, num_layers, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Nutrient limitation
        self.nutrient_limitation = tracer["parameters"]["nutrient_limitation"]
        self.nutrient_limitation_factor = {}
        self.nutrient_colimitation = 0.

        # Temperature regulation
        self.temperature_regulation = tracer["parameters"]["temperature_regulation"]
        self.temp_regulation_factor = 1.

        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Bacteria: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p','fe']
                if key in available_elements:
                    # Add constituent to composition/concentration
                    self.composition.append(key)

                    # Set initial conditions
                    if isinstance(tracer["composition"][key], str): # Read initial conditions from file
                        conc.append( np.fromfile(os.getcwd() + tracer["composition"][key]) )
                    elif isinstance(tracer["composition"][key], (int,float)): # Create array of initial conditions
                        if num_layers == 1: # 0d configuration
                            conc.append(tracer["composition"][key])
                        else: # 1d configuration
                            conc.append(tracer["composition"][key] * np.ones(num_layers))
                    elif isinstance(tracer["composition"][key], (list,np.ndarray)):
                        conc.append(np.array(tracer["composition"][key]))
                    elif isinstance(tracer["composition"][key], dict): # Create array of initial conditions based off ratio to base element
                        index = self.composition.index(base_element)
                        if num_layers == 1: # 0d configuration
                            conc.append(tracer["composition"][key][base_element] * conc[index] )
                        else: # 1d configuration
                            conc.append(tracer["composition"][key][base_element] * conc[index] * np.ones(num_layers))

                else:
                    sys.exit("Bacteria: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        # hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        # hold[...,0] = conc
        # self.conc = hold
        # self.d_dt = np.zeros_like(conc)
        # self.conc_ratio = np.zeros_like(conc)

        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.conc = np.zeros((len(self.composition),num_layers-1,iters),dtype=float)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const][:-1]
        else:   # Model as single box
            self.conc = np.zeros((len(self.composition),iters),dtype=float)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const]
        self.d_dt = np.zeros_like(self.conc[...,0],dtype=float)
        self.conc_ratio = np.ones_like(self.conc[...,0],dtype=float)
        

        # Production
        self.upt = {}   # Uptake

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
        
        # Reorder reactions (uptake needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]


    def bac(self, iter, base_element, physical, tracers):
        
        # Calculate oxygen limitation factor (if necessary)
        if "oxygen_inhibition" in self.__dict__:
            # Optional minimum concentration for aerobic/anerobic operations
            if "min_o2" in self.oxygen_inhibition:  o2 = np.maximum(tracers["o2"].conc[...,iter] , self.oxygen_inhibition["min_o2"] * np.ones_like(tracers["o2"].conc[...,iter]))
            else:   o2 = tracers["o2"].conc[...,iter]
            self.oxy_limitation_factor = monod(o2, self.oxygen_inhibition["half_sat"], self.oxygen_inhibition["exponent"])

        pass


    def add_nutrient(self, nutrient):
        """
        Add nutrients to phytoplankton and append dictionary of uptake rates
        """
        self.upt[nutrient] = np.zeros_like(self.conc[0,...],dtype=float)
        self.nutrient_limitation_factor[nutrient] = np.zeros_like(self.conc[0,...],dtype=float)

    
    def calculate_nutrient_limitation(self, base_element, iter, tracers):
        """
        Definition:: Calculates nutrient limitation factor as either a minimum, product, or sum of all nutrients which limit phytoplankton growth
        """

        fN = []
        for key in self.nutrient_limitation:
            if key != "colimitation":
                # Get nutrient chemical constituent
                element = tracers[key].composition[0]

                # Get index of nutrient chemcical constituent in phytoplankton composition dictionary
                element_index = self.composition.index(element)

                # Get index of nutrient chemcical constituent in cell quota dictionary
                quota_index = self.cell_quota.index(element)
                
                if self.nutrient_limitation[key]["type"] == "internal":
                    # Calculate nutrient limitation factor
                    func = ( self.conc_ratio[element_index] - self.cell_quota["min"][quota_index]) / ( self.cell_quota["opt"][quota_index] - self.cell_quota["min"][quota_index] )
                    
                    # Ensures nonzero value
                    func = np.maximum(1.E-20*np.ones_like(func), func)

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = func

                    # Append fN for colimitation calculation
                    if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                        fN.append(func)
                
                elif self.nutrient_limitation[key]["type"] == "external":
                    # Determine if Hill exponent exists for monod function
                    if "exponent" in self.nutrient_limitation[key]:
                        exponent = self.nutrient_limitation[key]["exponent"]
                    else:   # Default to 1. (no scaling)
                        exponent = 1.

                    # Calculate nutrient limitation factor
                    func = monod(tracers[key].conc[...,iter], self.nutrient_limitation[key]["half_sat"], exponent)
                    
                    # Ensures nonzero value
                    func = np.maximum(1.E-20*np.ones_like(func), func)

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = func

                    # Append fN for colimitation calculation
                    if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                        fN.append(func)

                else:
                    sys.exit("Nutrient limitation type not recognized. Check documentation and edit input file.")

        if "colimitation" in self.nutrient_limitation.keys():   # Multiple nutrients available
            fN = np.array(fN)
            if self.nutrient_limitation["colimitation"] == "minimum":
                self.nutrient_colimitation = np.min(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "product":
                self.nutrient_colimitation = np.prod(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "sum":
                self.nutrient_colimitation == np.sum(fN,  axis=0)
            else:
                sys.exit("Nutrient colimitation not recognized. Check documentation and edit input file.")
        else:   # Only one nutrient available
            self.nutrient_colimitation = fN

        
    def lysis(self, iter, base_element, parameters, c, p, ec, ep, ic, ip, tracers):

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Locate index of base element
        index = self.composition.index(base_element)
        bac = np.array(tracers[self.abbrev].conc[index][iter])

        lysis = parameters["lysis_rate"] * self.temp_regulation_factor * ( bac**2 )

        # Convert lysis rate (if necessary)
        if "convert_lysis" in parameters:
            if parameters["convert_lysis"] == "cell_quota":
                quota_index = self.nutrient_limitation["nutrients"].index(c)
                lysis *= self.nutrient_limitation_factor[quota_index]
            else:
                if isinstance(parameters["convert_lysis"],(int,float)) and not isinstance(parameters["convert_lysis"],bool):
                    lysis *= parameters["convert_lysis"]
                elif isinstance(parameters["convert_lysis"],str):
                    lysis *= float(Fraction(parameters["convert_lysis"]))
            
        
        # Update d_dt
        tracers[c].d_dt -= lysis

        # Apply partition to organic matter group (if necessary)
        if "partition" in parameters:   tracers[p].d_dt += lysis * parameters["partition"]
        else:                           tracers[p].d_dt += lysis


    def mortality(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates the non-grazing mortality of planktoninc species
        """

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        tc = np.array(tracers[c].conc[ic][iter])

        # Calculate mortality rate
        mortality = ( parameters["mortality_rate"][0] * tc ) + ( parameters["mortality_rate"][1] * (tc**2) )

        # Oxygen limitation
        if "oxygen_limited" in parameters and parameters["oxygen_limited"]:
            oxy_limitation_factor = np.minimum(1., nutrient_limitation(tracers["o2"].conc[...,iter], parameters["half_sat_oxygen"]))
            mortality += (1. - oxy_limitation_factor) * parameters["mortality_rate_oxy"] * tc

        # # Temperature regulation
        # if self.temp_limited:
        #     mortality = mortality #* self.temp_regulation_factor

        # Calculate concentration ratios
        # concentration_ratio(iter, ic, tracers[c])
        # concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        if parameters != None and "partition" in parameters:   # Mortality can be partitioned between dissolved and particulate detrital pools
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality * parameters["partition"]

        else:
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
            tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality


    def uptake(self, iter, base_element, parameters, c, p, ec, ep, ic, tracers):
        
        # Extract dict
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Identify the chemical constituent of the nutrient(s)
        element = ec[c[0]]

        # Get concentration of constituent in phytoplankton if present
        if element in self.composition:     element_index = self.composition.index(element)

        if parameters["strategy"] == "coupled":
            coupled_uptake = parameters["coupled_uptake"]
            linked_nutrients = coupled_uptake["link"]

            # Extract uptake rates of linked nutrients
            uptake_rates = []
            for nut in linked_nutrients:
                uptake_rates.append(self.uptake_rates[nut])

            # Calculate total linked uptake rate if multiple linked nutrients are used
            if len(linked_nutrients) > 1:   # use numpy "maximum" to ensure minimum uptake of 0.
                if coupled_uptake["method"] == "max":       linked_uptake = np.maximum(np.maximum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "min":     linked_uptake = np.maximum(np.minimum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "product": linked_uptake = np.maximum(np.prod(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "sum":     linked_uptake = np.maximum(np.sum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))

            if isinstance(parameters["convert_uptake"],(int,float)) and not isinstance(parameters["convert_uptake"],bool):
                uptake = self.nutrient_limitation_factor[element] * linked_uptake * coupled_uptake["convert_uptake"]
            elif isinstance(parameters["convert_uptake"],str):
                uptake = self.nutrient_limitation_factor[element] * linked_uptake * float(Fraction(parameters["convert_uptake"]))

            # Update d_dt
            tracers[c].d_dt -= np.array(ec) * np.maximum(uptake, np.zeros_like(uptake))
            if element in self.composition:     self.d_dt[element_index] += np.maximum(uptake, np.zeros_like(uptake))

        elif parameters["strategy"] == "independent":
            # Get concentration of element in bacterioplankton
            if element in self.composition:
                index = self.composition.index(element)
            else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                index = self.composition.index(base_element)

            bac = np.array(self.conc[index][iter])

            # Calculate maximum uptake rate
            max_uptake = (parameters["max_growth_rate"] + parameters["basal_metabolic_rate"]) / parameters["max_growth_efficiency"]

            # Calculate actual uptake
            uptake = max_uptake * monod(tracers[c].conc[ic][iter], self.nutrient_limitation_factor["dom1"], 1.) * bac

            # Multiply by conversion (if necessary)
            if "convert_uptake" in parameters:
                if isinstance(parameters["convert_uptake"],(int,float)) and not isinstance(parameters["convert_uptake"],bool):
                    uptake *= coupled_uptake["convert_uptake"]
                elif isinstance(parameters["convert_uptake"],str):
                    uptake *= float(Fraction(parameters["convert_uptake"]))
            
            if self.temperature_regulation["temp_limited"]:
                uptake *= self.temp_regulation_factor

            # Update d_dt
            tracers[c].d_dt -= uptake
            if element in self.composition:     self.d_dt[index] += uptake



    def respiration(self, iter, base_element, parameters, c, p, tracers):

        # Locate index of base element
        index = self.composition.index(base_element)

        # Get concentration of base element
        bac = tracers[self.abbrev].conc[index[iter]]

        # Calculate respiration rate
        respiration = self.temp_regulation_rator * parameters["respiration_rate"] * bac
        
        # Aeorbic --> O2 repired
        if tracers["o2"].conc[...,iter] > self.oxygen_inhibition["min_o2"]:
            # tracers["o2"].d_dt -= respiration * parameters["convert_o2"]
            if "o2" in c:   
                if isinstance(parameters["convert_o2"],(int,float)) and not isinstance(parameters["convert_o2"],bool):
                    tracers["o2"].d_dt -= respiration * parameters["convert_o2"]
                elif isinstance(parameters["convert_o2"],str):
                    tracers["o2"].d_dt -= respiration * float(Fraction(parameters["convert_o2"]))

        if "co2" in p:  
            if base_element == "c":     tracers["co2"].d_dt += respiration  
            else:                       
                if isinstance(parameters["convert_co2"],(int,float)) and not isinstance(parameters["convert_co2"],bool):
                    tracers["co2"].d_dt += respiration  * parameters["convert_co2"]   
                elif isinstance(parameters["convert_co2"],str):   
                    tracers["co2"].d_dt += respiration  * float(Fraction(parameters["convert_co2"]))

        # Update d_dt
        self.d_dt[index] -= respiration
        tracers[p].d_dt += respiration



















        