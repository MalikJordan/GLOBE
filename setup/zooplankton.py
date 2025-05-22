import copy
import os
import sys
import numpy as np
from functions.other_functions import concentration_ratio, nutrient_limitation, temperature_regulation, tracer_elements

class Zooplankton():
    """
    """

    def __init__(self, abbrev, iters, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Intracellular nutrinet quotas
        if "cell_quota" in tracer["parameters"]:
            self.cell_quota = tracer["parameters"]["cell_quota"]
        
        # Temperature regulation
        self.temp_limited = tracer["parameters"]["temp_limited"]
        if self.temp_limited:
            self.q10 = tracer["parameters"]["q10"]
        self.temp_regulation_factor = 1.

        # Oxygen limitation
        self.oxygen_limited = tracer["parameters"]["oxygen_limited"]
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
                    self.composition.append(key)
                    # conc.append(np.array(tracer["composition"][key]))
                    conc.append(tracer["composition"][key])
                else:
                    sys.exit("Zooplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
        hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        for i in range(0,len(conc)):
            hold[i,0] = np.array(conc[i])
        # hold[...,0] = conc
        self.conc = np.array(hold)
        self.d_dt = np.zeros_like(conc)
        # self.conc_ratio = np.zeros_like(conc)
        self.conc_ratio = np.copy(self.cell_quota["opt"])

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


    def zoo(self, iter, base_element, base_temp, temperature, tracers):
        
        check_conc = self.conc[:,iter]
        # Zero out grazing rates
        self.grazing_rates = {prey: 0.
                              for prey in self.grazing_rates}
        
        # Calculate temp regulation factor        
        if self.temp_limited:
            self.temp_regulation_factor = temperature_regulation(base_temp, temperature, self.q10)

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
                self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basl_respiration)
            if reac["type"] == "mortality":     self.mortality(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "respiration":   activity_respiration, basl_respiration = self.respiration(iter, reac["parameters"], c, p, tracers)
        
        if iter % 50 == 0:
            x=1
        x=1
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


    def excretion(self, iter, parameters, c, p, ec, ep, ic, ip, tracers, all_grazing, activity_respiration, basal_respiration):
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
        nutrient_index = list(ec).index(1.0)
        tc = np.array(tracers[c].conc[nutrient_index][iter])
        tp = np.array(tracers[p].conc[ip][iter])

        if "c" in self.composition:
            carbon_index = self.composition.index("c")

        if parameters["function"] == "constant":
            excretion = parameters["excretion_rate"] * tc

            # Calculate excretion of excess nutrient (above optimal nutrient quota) if necessary
            if tracers[p].type == "inorganic":
                if "c" in self.composition:
                    carbon_index = self.composition.index("c")
                    carbon_ratio = tc / np.array(tracers[c].conc[carbon_index][iter])

                    excretion = excretion * np.maximum(0,carbon_ratio - parameters["optimal_nutrient_quota"])
            
        elif parameters["function"] == "grazing":
            # Sum grazing rates for all chemical constituents
            graze_sum = np.zeros(len(self.composition))
            for const in self.composition:
                const_index = self.composition.index(const)
                graze_sum[const_index] = np.sum(all_grazing[const_index])
            
            # # Sum grazing rates for base element
            # graze_sum = sum([rate[ip] # * tracers[prey].conc[ic][iter]
            #                 for prey, rate in self.grazing_rates.items()])

            # excretion = self.assimilation_efficiency * ( 1 - self.ingestion_efficiency ) * graze_sum
            if tracers[p].type == "detritus":
                # excretion = self.ingestion_efficiency * ( 1 - self.assimilation_efficiency ) * graze_sum
                excretion = self.ingestion_efficiency * graze_sum
                if "c" in self.composition: # Scale excretion of carbon by assimilation efficiency
                    excretion[carbon_index] = excretion[carbon_index] * ( 1. - self.assimilation_efficiency )

            elif tracers[p].type == "inorganic":
                excreted_carbon = np.maximum(np.zeros_like(graze_sum[nutrient_index]), graze_sum[carbon_index] * (1. - self.ingestion_efficiency) - activity_respiration)
                excreted_nutrient = np.maximum(np.zeros_like(graze_sum[nutrient_index]), ( all_grazing[nutrient_index] * (1. - self.ingestion_efficiency) ) + ( basal_respiration * self.conc_ratio[nutrient_index] ))
                # graze_nutrient = (parameters["optimal_nutrient_quota"] * graze_sum)/(graze_sum + 1.E-20) - parameters["optimal_nutrient_quota"]
            
                excretion = np.maximum(np.zeros_like(excreted_carbon), excreted_nutrient/(excreted_carbon + 1.E-20) - self.cell_quota["opt"][nutrient_index] ) * excreted_carbon

        # Update d_dt
        if tracers[p].type == "detritus": # Scale by concentration ratio for excretion to detrital pools
            if parameters != None and "partition" in parameters:
                all_excretion = ec * excretion * parameters["partition"]
                tracers[c].d_dt -= ec * excretion * parameters["partition"]
                tracers[p].d_dt += ep * excretion * parameters["partition"]
            else:
                tracers[c].d_dt -= ec * excretion
                tracers[p].d_dt += ep * excretion

        else:
            tracers[c].d_dt -= ec * excretion
            tracers[p].d_dt += np.array(ep) * excretion
    

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
        tc = np.array(tracers[c].conc[ic][iter])
        tp = np.array(tracers[p].conc[ip][iter])

        # Extract grazing preference for prey
        pref = self.grazing_preferences[c] * tc

        # Calculate half saturation constant for grazing
        if parameters["function"] != "ivlev":
            # Calculate total food availability
            for prey, preference in self.grazing_preferences.items():
                # Concentration of base element in prey
                conc_prey = tracers[prey].conc[ic][iter]

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
            
            grazing = ( parameters["max_grazing_rate"] * self.prey_availability[c] ) / ( total_available + parameters["half_sat_grazing"] ) * tp
            
        elif parameters["function"] == "holling-3": # Sigmoidal
            # Calculate grazing preference
            pref = self.grazing_preferences[c] * ( tc**2 )
            pref_sum = sum([preference * ( tracers[prey].conc[ic][iter]**2 )
                            for prey, preference in self.grazing_preferences.items()])
            grazing = ( parameters["max_grazing_rate"] * pref ) / ( ( parameters["half_sat_grazing"]**2 ) + pref_sum ) * tp

            # grazing = ( parameters["max_grazing_rate"] * pref ) / ( ( half_sat_grazing**2 ) + pref_sum ) * tp
        
        # Temperature regulation
        if self.temp_limited:
            grazing = grazing * self.temp_regulation_factor

        rugc = self.temp_regulation_factor * parameters["max_grazing_rate"] * nutrient_limitation(self.prey_availability[c], parameters["half_sat_grazing"]) * tp
        sut = rugc / self.prey_availability[c]
        grazing = sut * self.prey_availability[c]


        gn = grazing * tracers[c].conc_ratio[1]
        gp = grazing * tracers[c].conc_ratio[2]
        gl = grazing * tracers[c].conc_ratio[3]

        phyto_conc_ratios = tracers[c].conc_ratio
        # Extract cell quotas from prey
        ratios = np.zeros(len(self.conc_ratio))
        for const in self.composition:
            if const in tracers[c].composition:
                index_prey = tracers[c].composition.index(const)
                index_pred = self.composition.index(const)
                # ratios[index_pred] = tracers[c].cell_quota["opt"][index_prey]
                ratios[index_pred] = tracers[c].conc_ratio[index_prey]

        all_grazing = ep * ratios * grazing   
        grazing_on_prey = ec * tracers[c].conc_ratio * grazing
        # grazing_on_phyto = ec * tracers[c].cell_quota["opt"] * grazing

        # Update d_dt
        tracers[c].d_dt -= grazing_on_prey
        tracers[p].d_dt += all_grazing
    
        # Update grazing rate dictionary
        tracers[p].grazing_rates[c] = all_grazing

        return all_grazing


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
            # tracers[c].d_dt -= ec * tracers[c].cell_quota["opt"] * mortality * parameters["partition"]
            # tracers[p].d_dt += ep * tracers[c].cell_quota["opt"] * mortality * parameters["partition"]
            all_mortality = ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality * parameters["partition"]

        else:
            # tracers[c].d_dt -= ec * tracers[c].cell_quota["opt"] * mortality
            # tracers[p].d_dt += ep * tracers[c].cell_quota["opt"] * mortality

            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
            tracers[p].d_dt += ep * tracers[c].conc_ratio * mortality


    def respiration(self, iter, parameters, c, p, tracers):
        """
        Definition:: Calculates zooplankton respiration
        """
        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        # zoo = self.conc[carbon_index][iter]
        zoo = tracers[self.abbrev].conc[carbon_index][iter]

        # Sum grazing rates for carbon constituent
        graze_sum = 0
        for prey, rate in self.grazing_rates.items():
            if "c" in tracers[prey].composition:    graze_sum += rate[tracers[prey].composition.index("c")]

        activity_respiration = (1 - self.assimilation_efficiency - self.ingestion_efficiency) * graze_sum
        basl_respiration = self.temp_regulation_factor * parameters["respiration_rate"] * zoo
        
        # # Temperature regulation
        # if self.temp_limited:   respiration += parameters["respiration_rate"] * self.temp_regulation_factor * zoo
        # else:   respiration += parameters["respiration_rate"] * zoo
        # basl_respiration = self.temp_regulation_factor * parameters["respiration_rate"] * zoo
        # Update d_dt
        total_respiration = activity_respiration + basl_respiration
        self.d_dt[carbon_index] -= total_respiration
        if "o2" in c:   tracers["o2"].d_dt -= total_respiration / parameters["mw_carbon"]
        if "co2" in p:  tracers["co2"].d_dt += total_respiration

        return activity_respiration, basl_respiration