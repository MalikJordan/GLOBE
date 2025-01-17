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
        
        # Temperature regulation
        self.temp_limited = tracer["parameters"]["temp_limited"]
        if self.temp_limited:
            self.q10 = tracer["parameters"]["q10"]
        self.fT = 1.

        # Grazing parameters
        self.grazing_preferences = tracer["parameters"]["grazing_preferences"]
        self.assimilation_efficiency = tracer["parameters"]["assimilation_efficiency"]
        self.ingestion_efficiency = tracer["parameters"]["ingestion_efficiency"]
        self.grazing_rates = {}
        
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
        self.conc_ratio = np.zeros_like(conc)

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
            
        # Reorder reactions (grazing needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "grazing"] + [item for item in self.reactions if item["type"] != "grazing"]


    def zoo(self, iter, base_element, base_temp, temperature, tracers):
        # Zero out grazing rates
        self.grazing_rates = {prey: 0.
                              for prey in self.grazing_rates}
        
        # Calculate temp regulation factor        
        if self.temp_limited:
            self.fT = temperature_regulation(base_temp, temperature, self.q10)

        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "grazing":       self.grazing(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "egestion":      self.egestion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "excretion":     self.excretion(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "mortality":     self.mortality(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "respiration":   rsp = self.respiration(iter, reac["parameters"], c, p, tracers)

        # return rsp
    

    def add_prey(self, prey):
        self.grazing_rates[prey] = 0.


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
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        if iter % 10000 == 0:
            q=1

        # Update d_dt
        if parameters != None and "partition" in parameters:   # Egestion can be partitioned between dissolved and particulate detrital pools
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * egestion * parameters["partition"]
            tracers[p].d_dt += ep * tracers[p].conc_ratio * egestion * parameters["partition"]

        else:
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * egestion
            tracers[p].d_dt += ep * tracers[p].conc_ratio * egestion


    def excretion(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates excretion of zooplankton to inorganic nutrient pool.
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

        if parameters["function"] == "constant":
            excretion = parameters["excretion_rate"] * tc

            # Calculate excretion of excess nutrient (above optimal nutrient quota) if necessary
            if "c" in self.composition:
                carbon_index = self.composition.index("c")
                carbon_ratio = tc / np.array(tracers[c].conc[carbon_index][iter])

                if tc == 0: 
                    q=1
                excretion = excretion * np.maximum(0,carbon_ratio - parameters["optimal_nutrient_quota"])
            
        elif parameters["function"] == "grazing":
            # Sum grazing rates for base element
            graze_sum = sum([rate[ip] # * tracers[prey].conc[ic][iter]
                            for prey, rate in self.grazing_rates.items()])

            excretion = self.assimilation_efficiency * ( 1 - self.ingestion_efficiency ) * graze_sum

        # Calculate concentration ratios
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        tracers[c].d_dt -= ec * tracers[c].conc_ratio * excretion
        tracers[p].d_dt += ep * tracers[p].conc_ratio * excretion
    

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

        # Calculate grazing function
        if parameters["function"] == "ivlev":       # Exponential
            function = parameters["max_grazing_rate"] * ( 1 - np.exp( -parameters["ivlev"] * tc ) ) * tp

        elif parameters["function"] == "holling-1": # Linear
            # Calculate slope based on prey concentration 
            if tc.size > 1:
                slope = np.zeros_like(tc)
                for i in range(0,len(slope)):
                    if tc[i] < (2 * parameters["half_sat_grazing"]):    slope[i] = parameters["max_grazing_rate"]/(2 * parameters["half_sat_grazing"])
                    else:   slope[i] = parameters["max_grazing_rate"]
            else:
                if tc < (2 * parameters["half_sat_grazing"]):    slope = parameters["max_grazing_rate"]/(2 * parameters["half_sat_grazing"])
                else:   slope = parameters["max_grazing_rate"]
            
            function = slope * tc * tp

        elif parameters["function"] == "holling-2": # Hyperbolic
            # Calculate grazing preference
            pref = self.grazing_preferences[c] * tc
            pref_sum = sum([preference * tracers[prey].conc[ic][iter]
                            for prey, preference in self.grazing_preferences.items()])

            function = ( parameters["max_grazing_rate"] * pref ) / ( parameters["half_sat_grazing"] + pref_sum ) * tp

        elif parameters["function"] == "holling-3": # Sigmoidal
            # Calculate grazing preference
            pref = self.grazing_preferences[c] * ( tc**2 )
            pref_sum = sum([preference * ( tracers[prey].conc[ic][iter]**2 )
                            for prey, preference in self.grazing_preferences.items()])

            function = ( parameters["max_grazing_rate"] * pref ) / ( ( parameters["half_sat_grazing"]**2 ) + pref_sum ) * tp
        
        # Temperature regulation
        if self.temp_limited:
            function = function * self.fT
        
        grazing_prey = function
        # grazing_predator = self.assimilation_efficiency * self.ingestion_efficiency * function
        grazing_predator = function

        # Calculate concentration ratios
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        tracers[c].d_dt -= ec * tracers[c].conc_ratio * grazing_prey
        # tracers[p].d_dt += ep * tracers[p].conc_ratio * grazing_predator
        tracers[p].d_dt += ec * tracers[p].conc_ratio * grazing_predator   # fix this to be generic if it works
        # tracers[p].d_dt += ec[:-1] * tracers[p].conc_ratio * grazing_predator   # fix this to be generic if it works
    

        # Update grazing rate dictionary
        tracers[p].grazing_rates[c] = grazing_prey * tracers[p].conc_ratio


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

        # Temperature regulation
        if self.temp_limited:
            mortality = mortality * self.fT

        # Oxygen limitation
        if "oxygen_limited" in parameters and parameters["oxygen_limited"]:
            # fO = nutrient_limitation(tc,parameters["half_sat_oxygen"])
            fO = tc / (parameters["half_sat_oxygen"] + tc)
            mortality = mortality * (1 - fO)

        # Calculate concentration ratios
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        if parameters != None and "partition" in parameters:   # Mortality can be partitioned between dissolved and particulate detrital pools
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality * parameters["partition"]

        else:
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
            tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality


    def respiration(self, iter, parameters, c, p, tracers):
        """
        Definition:: Calculates zooplankton respiration
        """
        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        zoo = self.conc[carbon_index][iter]

        # Sum grazing rates for carbon constituent
        graze_sum = 0
        for prey, rate in self.grazing_rates.items():
            if "c" in tracers[prey].composition:    graze_sum += rate[tracers[prey].composition.index("c")]

        respiration = (1 - self.assimilation_efficiency - self.ingestion_efficiency) * graze_sum
        
        # Temperature regulation
        if self.temp_limited:   respiration += parameters["respiration_rate"] * self.fT * zoo
        else:   respiration += parameters["respiration_rate"] * zoo

        # Update d_dt
        self.d_dt[carbon_index] -= respiration
        if "o2" in c:   tracers["o2"].d_dt -= respiration / parameters["mw_carbon"]
        if "co2" in p:  tracers["co2"].d_dt += respiration

        return respiration