import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, nutrient_limitation, temperature_regulation, tracer_elements

class Inorganic():
    """
    
    """

    def __init__(self, long_name, composition, type):
        self.name = long_name
        self.conc = []
        self.type = type

        if len(composition) > 1:
            sys.exit("Inorganic: Only one element accepted per inorganic nutrient. Check documentation adn edit input file.")
        elif len(composition) < 1:
            sys.exit("Inorganic: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition: 
                self.conc.append(composition[key])
            self.d_dt = np.zeros_like(self.conc)

    def nitrification(self, parameters, fT, fO, tracers):
        nitrification = fT * fO * parameters["nitrification_rate"] * tracers["nh4"].conc

        tracers["nh4"].d_dt -= nitrification
        if "o2" in tracers: tracers["o2"].d_dt -= 2. * nitrification    # '2 *' for stoichiometry, NH4(+) + 2O2 --> NO3(-) + H2O + 2H(+)
        tracers["no3"].d_dt += nitrification
        
    
    def inorg(self, base_element, base_temp, coordinates, mixed_layer_depth, reactions, surface_PAR, temperature, tracers):
        
        for reac in reactions:
            tc = tracers[reac["consumed"]]  # consumed tracer
            tp = tracers[reac["produced"]]  # produced tracer

            # ic/ip = index of base element in consumed/produced tracer
            # ec/ep = array of  elements in consumed/produced tracer affected by current reaction
            ec, ep, ic, ip = tracer_elements(base_element, reac, tc, tp)


            if reac["type"] == "nitrification":  
                fT = temperature_regulation(base_temp, temperature, reac["parameters"]["q10"])   # calculate temperature regulation factor for nitrification
                if "o2" in tracers: fO = nutrient_limitation(tracers["o2"].conc,tracers["o2"].nutrient_half_sat)
                else:   fO = 1.
                self.nitrification(reac["parameters"], tracers)
            