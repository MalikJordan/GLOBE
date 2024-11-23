import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Detritus():
    """
    
    """

    def __init__(self, abbrev, iters, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]
        
        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Detritus: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    self.composition.append(key)
                    conc.append(tracer["composition"][key])
                else:
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
        hold = np.zeros((len(conc),iters))
        hold[...,0] = conc
        self.conc = np.array(hold)
        self.d_dt = np.zeros_like(conc)
        self.conc_ratio = np.zeros_like(conc)
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
          

    def detritus(self, iter, base_element, tracers):
    
        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            
            if reac["type"] == "remineralization":  self.remineralization(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
    
    
    def remineralization(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        
        # Extract dict
        if len(c) > 1 and "o2" in c:
        # Remineralization of carbon also affects oxygen (sink)
            for t in c:
                if t == "o2":   pass
                else:
                    consumed = t
                    break
        else:   consumed = c[0]
        ec = ec[consumed]
        ic = ic[consumed]

        # Get concentration of remineralized nutrient in organic matter pool
        index = list(ec).index(1.)
        tc = np.array(tracers[consumed].conc[index][iter])

        if p[0] == None:    pass
        else:
            p = p[0]
            ep = ep[p]
            ip = ip[p]
            tp = np.array(tracers[p].conc[ip][iter])
        
        remineralization = (parameters["remineralization_rate"]) * tc

        # concentration_ratio(iter, ic, tracers[consumed])
        # tracers[consumed].d_dt -= ec * tracers[consumed].conc_ratio * remineralization
        tracers[consumed].d_dt -= ec * remineralization
        if "o2" in c:
            tracers["o2"].d_dt -=  remineralization / parameters["mw_carbon"]
        
        if p[0] == None:    pass
        else:   
            # concentration_ratio(iter, ip, tracers[p])
            # tracers[p].d_dt += ep * tracers[p].conc_ratio * remineralization
            tracers[p].d_dt += np.array(ep) * remineralization

    
    