import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Detritus():
    """
    
    """

    def __init__(self, long_name, composition):
        self.name = long_name
        self.composition = []
        self.conc = []
        self.type = "detritus"
        
        if len(composition) < 1:
            sys.exit("Detritus: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition:
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    self.composition.append(key)
                    self.conc.append(composition[key])
                else:
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
    
        self.d_dt = np.zeros_like(self.conc)
        self.conc_ratio = np.zeros_like(self.conc)

    
    def remineralization(self, parameters, ec, ep, ic, ip, consumed, produced):

        remineralization = (parameters["remineralization_rate"]) * np.array(consumed.conc[ic])

        concentration_ratio(ic,consumed)
        concentration_ratio(ip,produced)

        consumed.d_dt[ic] -= ec * consumed.conc_ratio * remineralization
        produced.d_dt[ip] += ep * produced.conc_ratio * remineralization

    
    def detritus(self, base_element, reactions, tracers):

        # Calculate bgc rates
        for reac in reactions:
            tc = tracers[reac["consumed"]]  # consumed tracer
            tp = tracers[reac["produced"]]  # produced tracer

            # ic/ip = index of base element in consumed/produced tracer
            # ec/ep = array of  elements in consumed/produced tracer affected by current reaction
            ec, ep, ic, ip = tracer_elements(base_element, reac, tc, tp)

            if reac["type"] == "remineralization":  self.remineralization(reac["parameters"], ec, ep, ic, ip, tc, tp)
