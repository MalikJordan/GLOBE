import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, nutrient_limitation, temperature_regulation, tracer_elements

class Inorganic():
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
        self.temp_regulation_factor = 1.

        # Oxygen regulation
        

        # Concentration array
        self.composition = []
        conc = []
        if len(tracer["composition"]) > 1:    sys.exit("Inorganic: Only one element accepted per inorganic nutrient. Check documentation adn edit input file.")
        elif len(tracer["composition"]) < 1:  sys.exit("Inorganic: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:   pass

        for key in tracer["composition"]: 
            self.composition.append(key)
            conc.append(tracer["composition"][key])
        hold = np.zeros((len(conc),iters))
        hold[...,0] = conc
        self.conc = np.array(hold)
        self.d_dt = np.zeros_like(conc)
        self.conc_ratio = np.ones_like(self.conc[...,0])

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
        
    
    def inorg(self, iter, base_element, base_temp, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, salinity, wind, tracers):
        check_conc = self.conc[:,iter]
        
        if self.temp_limited:
            self.temp_regulation_factor = temperature_regulation(base_temp, temperature, self.q10)   # calculate temperature regulation factor for nitrification
        
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "nitrification" and self.abbrev == "no3":  
                if "o2" in tracers: oxy_limitation_factor = nutrient_limitation(tracers["o2"].conc[...,iter],reac["parameters"]["half_sat_oxygen"])
                else:   oxy_limitation_factor = 1.
                self.nitrification(iter, reac["parameters"], oxy_limitation_factor, tracers)
            if reac["type"] == "reaeration":    self.reaeration(iter, reac["parameters"], dz, temperature, salinity, wind)
            if reac["type"] == "reoxidation" and self.abbrev == "hs":   
                if "o2" in tracers: oxy_limitation_factor = nutrient_limitation(tracers["o2"].conc[...,iter],reac["parameters"]["half_sat_oxygen"])
                else:   oxy_limitation_factor = 1.
                self.reoxidation(iter, reac["parameters"], oxy_limitation_factor, tracers)
        
        if iter % 50 == 0:
            x=1
        if self.abbrev == 'nh4':
            x=1
        if self.abbrev == 'o2':
            x=1
        x=1

    def nitrification(self, iter, parameters, oxy_limitation_factor,  tracers):
        """
        Definition:: Calculates nitrification rate
        """
        nitrification = self.temp_regulation_factor * oxy_limitation_factor * parameters["nitrification_rate"] * tracers["nh4"].conc[...,iter]
        nitrification = np.maximum(np.zeros_like(nitrification),nitrification)
        
        tracers["nh4"].d_dt -= nitrification
        if "o2" in tracers: tracers["o2"].d_dt -= parameters["nitrification_stoic_coeff"] * nitrification    # '2. *' for stoichiometry, NH4(+) + 2O2 --> NO3(-) + H2O + 2H(+)
        tracers["no3"].d_dt += nitrification


    def reaeration(self, iter, parameters, dz, temperature, salinity, wind):

        # Calculate absolute temperature divided by 100
        abt = (temperature + 273.15) / 100.

        # Calculate theoretical oxygen saturatino for temp + salt and conver into proper units of [mmol O2 / m3]
        # oxy_sat = np.exp(-173.4292 + (249.6339/abt) + (143.3483*np.log(abt))-(21.8492*abt) + salinity*(-0.033096 + 0.014259*abt - 0.0017*(abt**2)))/(24.4665E-3)
        oxy_sat = np.exp(-173.4292 + (249.6339/abt) + (143.3483*np.log(abt))-(21.8492*abt) + salinity*(-0.033096 + 0.014259*abt - 0.0017*(abt**2)))*44.661

        # Calculate Schmidt number, ratio between the kinematic viscosity and the molecular diffusivity of CO2
        schmidt_number = parameters["k1"] - ( parameters["k2"]*temperature ) + ( parameters["k3"]*(temperature**2) ) - ( parameters["k4"]*(temperature**3) )
        schmidt_ratio = parameters["schmidt"] / schmidt_number

        # Schmidt ratio limited to 0
        schmidt_ratio = np.maximum(np.zeros_like(schmidt_ratio),schmidt_ratio)

        # Calculate wind dependency, includif conversion from cm/hr to m/s
        wind_dependency = ( parameters["d"]*(wind**2))*np.sqrt(schmidt_ratio)

        # Convert from cm/hr to m/day
        wind_dependency = wind_dependency * 0.01 * 24

        # Calculate flux of o2
        d_o2 = wind_dependency * (oxy_sat - self.conc[...,iter]) / dz

        # Update d_dt
        self.d_dt += d_o2


    def reoxidation(self, iter, parameters, oxy_limitation_factor, tracers):
        """
        Definition:: Calculation reoxidation rate of reduction equivalents
        """
        reoxidation = parameters["reoxidation_rate"] * oxy_limitation_factor * tracers["hs"].conc[..., iter]

        # Update d_dt
        tracers["hs"].d_dt -= reoxidation
        tracers["o2"].d_dt -= reoxidation / parameters["reoxidation_stoic_coeff"]