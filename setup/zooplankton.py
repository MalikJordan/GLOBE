import copy
import os
import sys
import numpy as np
from functions.other_functions import concentration_ratio, nutrient_limitation, temperature_dependence, tracer_elements
from fractions import Fraction
class Zooplankton():
    """
    """

    def __init__(self, abbrev, base_element, iters, num_layers, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Intracellular nutrinet quotas
        if "cell_quota" in tracer["parameters"]:
            self.cell_quota = tracer["parameters"]["cell_quota"]
        
        # Temperature regulation
        self.temperature_regulation = tracer["parameters"]["temperature_regulation"]
        self.temp_regulation_factor = 1.

        # Oxygen inhibition
        self.oxygen_inhibition = tracer["parameters"]["oxygen_inhibition"]
        self.oxy_limitation_factor = 1.

        # Grazing parameters
        self.grazing_preferences = tracer["parameters"]["grazing_preferences"]
        self.assimilation_efficiency = tracer["parameters"]["assimilation_efficiency"]
        self.ingestion_efficiency = tracer["parameters"]["ingestion_efficiency"]
        self.grazing_rates = {}
        self.prey_availability = {}
        
        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Zooplankton: Element required for " + self.name + ". Check documentation adn edit input file.")
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
                    sys.exit("Zooplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        # hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        # hold[...,0] = conc
        # self.conc = hold
        # self.d_dt = np.zeros_like(conc)
        # self.conc_ratio = np.copy(self.cell_quota["opt"])

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

        # Reorder reactions (grazing needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "respiration"] + [item for item in self.reactions if item["type"] != "respiration"]
        self.reactions = [item for item in self.reactions if item["type"] == "grazing"] + [item for item in self.reactions if item["type"] != "grazing"]


    def zoo(self, iter, base_element, physical, tracers):
        
        # check_conc = self.conc[:,iter]
        # Zero out grazing rates
        self.grazing_rates = {prey: 0.
                              for prey in self.grazing_rates}
        self.graze_sum = np.zeros_like(self.conc_ratio, dtype=float)
        
        # Calculate temp regulation factor (if necessary)
        if self.temperature_regulation["temp_limited"]:
            # self.temp_regulation_factor = temperature_dependence(base_temp, temperature, self)
            self.temp_regulation_factor = temperature_dependence(physical["bgc_phys_vars"]["temperature"], self)

        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "grazing":       all_grazing = self.grazing(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "egestion":      self.egestion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "excretion":     
                # if self.calc_respiration:   self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basl_respiration)
                # else:                       self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, 0., 0.)
                if self.calc_respiration:   pass
                else:
                    activity_respiration = 0.
                    basl_respiration = 0.
                if self.calc_grazing:   pass
                else:   all_grazing = np.zeros_like(self.conc_ratio)
                self.excretion(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basl_respiration)
            if reac["type"] == "mortality":     self.mortality(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "respiration":   activity_respiration, basl_respiration = self.respiration(iter, base_element, reac["parameters"], c, p, tracers)
        
        # if iter % 50 == 0:
        #     x=1
        # x=1
        # return rsp
    

    def add_prey(self, prey):
        self.grazing_rates[prey] = 0.
        self.prey_availability[prey] = 0.


    def egestion(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates zooplankton loss to detrital pool as feacal pellet production
        Return:: Egestion rate
        """

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Sum grazing rates for base element
        graze_sum = sum([rate[ip] #* tracers[prey].conc[ic][iter]
                        for prey, rate in self.grazing_rates.items()])

        # Calculate egestion
        egestion = ( 1 - self.assimilation_efficiency ) * graze_sum

        # Calculate concentration ratios
        # concentration_ratio(iter, ic, tracers[c])
        # concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        if parameters != None and "partition" in parameters:   # Egestion can be partitioned between dissolved and particulate detrital pools
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * egestion * parameters["partition"]
            tracers[p].d_dt += ep * tracers[p].conc_ratio * egestion * parameters["partition"]

        else:
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * egestion
            tracers[p].d_dt += ep * tracers[p].conc_ratio * egestion


    def excretion(self, iter, base_element, parameters, c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basal_respiration):
        """
        Definition:: Calculates excretion of zooplankton to nutrient or detrital pool.
                     Excretion can be represented either as a constant rate or as a fraction of zooplankton grazing on phytoplankton.
        Return:: Excretion rate
        """

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        nutrient_index = list(ec).index(1.0)    # Index of excreted nutrient
        tc = np.array(tracers[c].conc[nutrient_index,:,iter])
        index = self.composition.index(base_element)    # Index of base element

        if parameters["function"] == "constant":
            excretion = parameters["excretion_rate"] * tc

            # Calculate excretion of excess nutrient (above optimal nutrient quota) if necessary
            if tracers[p].type == "inorganic":
                element_ratio = tc / np.array(tracers[c].conc[index][iter])
                excretion = excretion * np.maximum(0,element_ratio - parameters["optimal_nutrient_quota"])

        elif parameters["function"] == "grazing":
            # Sum grazing rates for all chemical constituents
            # graze_sum = np.zeros(len(self.composition))
            # for const in self.composition:
            #     const_index = self.composition.index(const)
            #     # graze_sum[const_index] = np.sum(all_grazing[const_index])
            #     graze_sum[const_index] = np.sum(self.grazing_rates[const_index])
            
            if tracers[p].type == "detritus":
                # excretion = self.ingestion_efficiency * graze_sum
                excretion = self.ingestion_efficiency * self.graze_sum
                excretion[index] = excretion[index] * ( 1. - self.assimilation_efficiency )

            elif tracers[p].type == "inorganic":
                # excreted_base = np.maximum(np.zeros_like(graze_sum[nutrient_index]), graze_sum[index] * (1. - self.ingestion_efficiency) - activity_respiration)
                # excreted_nutrient = np.maximum(np.zeros_like(graze_sum[nutrient_index]), ( all_grazing[nutrient_index] * (1. - self.ingestion_efficiency) ) + ( basal_respiration * self.conc_ratio[nutrient_index] ))
                excreted_base = np.maximum(np.zeros_like(self.graze_sum[nutrient_index]), self.graze_sum[index] * (1. - self.ingestion_efficiency) - activity_respiration)
                # excreted_nutrient = np.maximum(np.zeros_like(self.graze_sum[nutrient_index]), ( all_grazing[nutrient_index] * (1. - self.ingestion_efficiency) ) + ( basal_respiration * self.conc_ratio[nutrient_index] ))
                excreted_nutrient = np.maximum(np.zeros_like(self.graze_sum[nutrient_index]), ( self.graze_sum[nutrient_index] * (1. - self.ingestion_efficiency) ) + ( basal_respiration * self.conc_ratio[nutrient_index] ))
            
                excretion = np.maximum(np.zeros_like(excreted_base), excreted_nutrient/(excreted_base + 1.E-20) - self.cell_quota["opt"][nutrient_index] ) * excreted_base
                # excretion = np.maximum(np.zeros_like(excreted_base), excreted_nutrient/(excreted_base + 1.E-20) - self.conc_ratio[nutrient_index] ) * excreted_base

        # Update d_dt
        if tracers[p].type == "detritus": # Scale by concentration ratio for excretion to detrital pools
            # Extract cell quotas from zooplankton
            # ratios = np.zeros(len(self.conc_ratio))
            ratios = np.zeros_like(tracers[p].conc_ratio)
            for const in self.composition:
                if const in tracers[p].composition:
                    index_zoo = tracers[c].composition.index(const)
                    index_om = self.composition.index(const)
                    ratios[index_om] = tracers[c].conc_ratio[index_zoo]

            if parameters != None and "partition" in parameters:
                # all_excretion = ec * excretion * parameters["partition"]
                # tracers[c].d_dt -= ec * excretion * parameters["partition"]
                # tracers[p].d_dt += ep * excretion * parameters["partition"]

                for i in range(len(ec)):
                    tracers[c].d_dt[i] -= ec[i] * excretion[i] * parameters["partition"][i] #* tracers[c].conc_ratio[i]
                for j in range(len(ep)):
                    tracers[p].d_dt[j] += ep[j] * excretion[j] * parameters["partition"][j] #* ratios[j]

            else:
                for i in range(len(ec)):
                    tracers[c].d_dt[i] -= ec[i] * excretion[i] * tracers[c].conc_ratio[i]
                for j in range(len(ep)):
                    tracers[p].d_dt[j] += ep[j] * excretion[j] * ratios[j]
                # tracers[c].d_dt -= ec * excretion
                # tracers[p].d_dt += ep * excretion

        else:
            for i in range(len(ec)):
                tracers[c].d_dt[i] -= ec[i] * excretion
            for j in range(len(ep)):
                tracers[p].d_dt[j] += ep[j] * excretion
            # tracers[c].d_dt -= ec * excretion
            # tracers[p].d_dt += np.array(ep) * excretion

    
    def grazing(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates the zooplankton grazing rate on a particular species using user choice of the Ivlev Equation,
                    Holling Type I Response, Holling Type II Response, or Holling Type III Response.
        Return:: grazing_c - Grazing rate to be added to "grazing_rates[]" for later summation
                grazing_p - Porition of grazing rate allocated to zooplaknton (scaled by assimilation and ingestion efficiencies).
        """

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        tc = np.array(tracers[c].conc[ic,:,iter])
        tp = np.array(tracers[p].conc[ip,:,iter])

        # Extract grazing preference for prey
        pref = self.grazing_preferences[c] * tc

        # Calculate half saturation constant for grazing
        if parameters["function"] != "ivlev":
            # Calculate total food availability
            for prey, preference in self.grazing_preferences.items():
                # Concentration of base element in prey
                conc_prey = tracers[prey].conc[ic,:,iter]

                # Capture efficiency for current prey in list of available
                eff_prey = conc_prey / ( conc_prey + parameters["feeding_threshold"] )

                # Total food availability is sum for all prey
                self.prey_availability[prey] = preference * eff_prey * conc_prey

        # Calculate grazing function
        if parameters["function"] == "ivlev":       # Exponential
            grazing = parameters["max_grazing_rate"] * ( 1 - np.exp( -parameters["ivlev"] * tc ) ) * tp

        elif parameters["function"] == "holling-1": # Linear
            # Calculate slope based on prey concentration 
            if tc.size > 1:
                slope = np.zeros_like(tc)
                for i in range(0,len(slope)):
                    if tc[i] < (2 * parameters["half_sat_grazing"]):    slope[i] = parameters["max_grazing_rate"]/(2 * parameters["half_sat_grazing"])
                    # if tc[i] < (2 * parameters["half_sat_grazing"]):    slope[i] = parameters["max_grazing_rate"]/(2 * half_sat_grazing)
                    else:   slope[i] = parameters["max_grazing_rate"]
            else:
                if tc < (2 * parameters["half_sat_grazing"]):    slope = parameters["max_grazing_rate"]/(2 * parameters["half_sat_grazing"])
                # if tc < (2 * parameters["half_sat_grazing"]):    slope = parameters["max_grazing_rate"]/(2 * half_sat_grazing)
                else:   slope = parameters["max_grazing_rate"]
            
            grazing = slope * tc * tp

        elif parameters["function"] == "holling-2": # Hyperbolic
            # Calculate total prey availability
            total_available = sum([avail for prey, avail in self.prey_availability.items()])
            
            # Calculate specific grazing rate for individual prey
            grazing = ( parameters["max_grazing_rate"] * self.prey_availability[c] ) / ( total_available + parameters["half_sat_grazing"] ) * tp

            # Calculate total uptake rate
            uptake = ( parameters["max_grazing_rate"] * total_available ) / ( total_available + parameters["half_sat_grazing"] ) * tp
            total_uptake = np.zeros_like(self.conc_ratio)
            for prey in self.prey_availability:
                for const in tracers[prey].composition:
                    if const in self.composition:
                        index_prey = tracers[prey].composition.index(const)
                        index_pred = self.composition.index(const)
                        total_uptake[index_pred] += (uptake / total_available) * tracers[prey].conc_ratio[index_prey] * self.prey_availability[prey]
                        # total_uptake[index_pred] += uptake * tracers[prey].conc_ratio[index_prey] * self.prey_availability[prey]
            
        elif parameters["function"] == "holling-3": # Sigmoidal
            # Calculate grazing preference
            pref = self.grazing_preferences[c] * ( tc**2 )
            pref_sum = sum([preference * ( tracers[prey].conc[ic][iter]**2 )
                            for prey, preference in self.grazing_preferences.items()])
            grazing = ( parameters["max_grazing_rate"] * pref ) / ( ( parameters["half_sat_grazing"]**2 ) + pref_sum ) * tp

            # grazing = ( parameters["max_grazing_rate"] * pref ) / ( ( half_sat_grazing**2 ) + pref_sum ) * tp
        
        # Temperature regulation
        if self.temperature_regulation["temp_limited"]:
            grazing = grazing * self.temp_regulation_factor
            total_uptake = total_uptake * self.temp_regulation_factor


        # Extract cell quotas from prey
        ratios = np.zeros_like(tracers[p].conc_ratio)
        for const in self.composition:
            if const in tracers[c].composition:
                index_prey = tracers[c].composition.index(const)
                index_pred = self.composition.index(const)
                ratios[index_pred] = tracers[c].conc_ratio[index_prey]

        # Update d_dt
        # grazing_on_prey = np.zeros_like(tracers[c].conc_ratio)
        for i in range(len(ec)):
            # grazing_on_prey[i] = ec[i] * tracers[c].conc_ratio[i] * grazing
            # tracers[c].d_dt[i] -= grazing_on_prey[i]
            tracers[c].d_dt[i] -= ec[i] * tracers[c].conc_ratio[i] * grazing

        all_grazing = np.zeros_like(tracers[p].conc_ratio)
        for j in range(len(ep)):
            all_grazing[j] = ep[j] * ratios[j] * grazing
            tracers[p].d_dt[j] += all_grazing[j]
  
        # Update grazing rate dictionary
        tracers[p].grazing_rates[c] = all_grazing
        
        # Update graze_sum
        for i in range(0,len(self.graze_sum)):
            # self.graze_sum[i] = self.conc_ratio[i] * total_uptake
            self.graze_sum[i] = total_uptake[i]
        # self.graze_sum += all_grazing

        sut = total_uptake / total_available
        x = tracers[c].conc_ratio * grazing
        # y = self.conc_ratio * total_uptake
        x = 1
    

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
        tc = np.array(tracers[c].conc[ic,:,iter])

        # Calculate mortality rate
        mortality = ( parameters["mortality_rate"][0] * tc ) + ( parameters["mortality_rate"][1] * (tc**2) )

        # Oxygen limitation
        if "oxygen_limited" in parameters and parameters["oxygen_limited"]:
            oxy_limitation_factor = np.minimum(1., nutrient_limitation(tracers["o2"].conc[...,iter], parameters["half_sat_oxygen"]))
            mortality += (1. - oxy_limitation_factor[0,:]) * parameters["mortality_rate_oxy"] * tc

        # # Temperature regulation
        # if self.temp_limited:
        #     mortality = mortality #* self.temp_regulation_factor

        # Calculate concentration ratios
        # concentration_ratio(iter, ic, tracers[c])
        # concentration_ratio(iter, ip, tracers[p])

        # Extract cell quotas from zooplankton
        # ratios = np.zeros(len(self.conc_ratio))
        ratios = np.zeros_like(tracers[p].conc_ratio)
        for const in self.composition:
            if const in tracers[p].composition:
                index_zoo = tracers[c].composition.index(const)
                index_om = self.composition.index(const)
                ratios[index_om] = tracers[c].conc_ratio[index_zoo]

        # Update d_dt
        if parameters != None and "partition" in parameters:   # Mortality can be partitioned between dissolved and particulate detrital pools
            # tracers[c].d_dt -= ec * tracers[c].cell_quota["opt"] * mortality * parameters["partition"]
            # tracers[p].d_dt += ep * tracers[c].cell_quota["opt"] * mortality * parameters["partition"]
            # all_mortality = ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            # tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            # tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality * parameters["partition"]

            # if parameters != None and "partition" in parameters:
                # all_excretion = ec * excretion * parameters["partition"]
                # tracers[c].d_dt -= ec * excretion * parameters["partition"]
                # tracers[p].d_dt += ep * excretion * parameters["partition"]

            for i in range(len(ec)):
                tracers[c].d_dt[i] -= ec[i] * mortality * parameters["partition"][i] * tracers[c].conc_ratio[i]
            for j in range(len(ep)):
                tracers[p].d_dt[j] += ep[j] * mortality * parameters["partition"][j] * ratios[j]

            # else:
            #     for i in range(len(ec)):
            #         tracers[c].d_dt[i] -= ec[i] * mortality
            #     for j in range(len(ep)):
            #         tracers[p].d_dt[j] += ep[j] * mortality
                # tracers[c].d_dt -= ec * excretion
                # tracers[p].d_dt += ep * excretion


        else:
            # tracers[c].d_dt -= ec * tracers[c].cell_quota["opt"] * mortality
            # tracers[p].d_dt += ep * tracers[c].cell_quota["opt"] * mortality

            # tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
            # tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality

            for i in range(len(ec)):
                tracers[c].d_dt[i] -= ec[i] * mortality * tracers[c].conc_ratio[i]
            for j in range(len(ep)):
                tracers[p].d_dt[j] += ep[j] * mortality * ratios[j]

        x=1
    

    def respiration(self, iter, base_element, parameters, c, p, tracers):
        """
        Definition:: Calculates zooplankton respiration
        """

        # Locate index of base element
        index = self.composition.index(base_element)

        # Get concentration of base element
        zoo = tracers[self.abbrev].conc[index,:,iter]

        # Sum grazing rates for base element
        # graze_sum = 0.
        # for prey, rate in self.grazing_rates.items():
        #     graze_sum += rate[tracers[prey].composition.index(base_element)]

        # activity_respiration = (1 - self.assimilation_efficiency - self.ingestion_efficiency) * graze_sum
        activity_respiration = (1 - self.assimilation_efficiency - self.ingestion_efficiency) * self.graze_sum[index]
        basal_respiration = self.temp_regulation_factor * parameters["respiration_rate"] * zoo

        # Update d_dt
        total_respiration = activity_respiration + basal_respiration
        self.d_dt[index] -= total_respiration
        if "o2" in c:   
            if isinstance(parameters["convert_o2"],(int,float)) and not isinstance(parameters["convert_o2"],bool):
                tracers["o2"].d_dt -= total_respiration * parameters["convert_o2"]
            elif isinstance(parameters["convert_o2"],str):
                tracers["o2"].d_dt -= total_respiration * float(Fraction(parameters["convert_o2"]))

        if "co2" in p:  
            if base_element == "c":     tracers["co2"].d_dt += total_respiration  
            else:                       
                if isinstance(parameters["convert_co2"],(int,float)) and not isinstance(parameters["convert_co2"],bool):
                    tracers["co2"].d_dt += total_respiration  * parameters["convert_co2"]   
                elif isinstance(parameters["convert_co2"],str):   
                    tracers["co2"].d_dt += total_respiration  * float(Fraction(parameters["convert_co2"]))

        return activity_respiration, basal_respiration
    