import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Detritus():
    """
    
    """

    def __init__(self, iters, reactions, **tracer):
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

        # Add relevant reactions
        self.reactions = []
        for reac in reactions:
            # Add reaction to dictionary
            if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
            else:   consumed = {"empty": "empty"}
            if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
            else:   produced = {"empty": "empty"}
            if ( list(tracer.keys())[0] in consumed.keys() ) or ( list(tracer.keys())[0] in produced.keys() ):
                self.reactions.append(reac)
          


    # def __init__(self, abbrev, composition, iters, long_name, reactions, type):
    #     self.name = long_name
    #     self.type = type
        
    #     # Composition and concentration arrays
    #     self.composition = []
    #     conc = []
    #     if len(composition) < 1:
    #         sys.exit("Detritus: Element required for " + long_name + ". Check documentation adn edit input file.")
    #     else:
    #         for key in composition:
    #             available_elements = ['c','n','p','chl','fe','si','caco3']
    #             if key in available_elements:
    #                 self.composition.append(key)
    #                 conc.append(composition[key])
    #             else:
    #                 sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
    #     hold = np.zeros((len(conc),iters))
    #     hold[...,0] = conc
    #     self.conc = np.array(hold)
    #     self.d_dt = np.zeros_like(conc)
    #     self.conc_ratio = np.zeros_like(conc)

    #     # Add relevant reactions
    #     self.reactions = []
    #     for reac in reactions:
    #         # Add reaction to dictionary
    #         if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
    #         else:   consumed = {"empty": "empty"}
    #         if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
    #         else:   produced = {"empty": "empty"}
    #         if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
    #             self.reactions.append(reac)
          

    
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
