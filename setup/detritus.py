import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements
from fractions import Fraction
class Detritus():
    """
    
    """


    def __init__(self, abbrev, base_element, iters, num_layers, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Light limitation
        if "light_attenuation" in tracer["parameters"]:
            self.light_attenuation = tracer["parameters"]["light_attenuation"]
        else:
            self.light_attenuation = 0.
        
        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Detritus: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p','chl','fe','si','caco3']
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
                    sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        hold[...,0] = conc
        self.conc = hold
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
          

    # def __init__(self, abbrev, iters, reactions, **tracer):
    #     self.abbrev = abbrev
    #     self.name = tracer["long_name"]
    #     self.type = tracer["type"]

    #     # Light limitation
    #     if "light_attenuation" in tracer["parameters"]:
    #         self.light_attenuation = tracer["parameters"]["light_attenuation"]
    #     else:
    #         self.light_attenuation = 0.
        
    #     # Composition and concentration arrays
    #     self.composition = []
    #     conc = []
    #     if len(tracer["composition"]) < 1:
    #         sys.exit("Detritus: Element required for " + self.name + ". Check documentation adn edit input file.")
    #     else:
    #         for key in tracer["composition"]:
    #             available_elements = ['c','n','p','chl','fe','si','caco3']
    #             if key in available_elements:
    #                 self.composition.append(key)
    #                 conc.append(tracer["composition"][key])
    #             else:
    #                 sys.exit("Detritus: Element '" + key + "' not recognized. Check documentation and edit input file.")
    #     hold = np.zeros((len(conc),iters))
    #     hold[...,0] = conc
    #     self.conc = np.array(hold)
    #     self.d_dt = np.zeros_like(conc)
    #     self.conc_ratio = np.zeros_like(conc)
    #     self.conc_ratio = np.ones_like(self.conc[...,0])

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
          

    def detritus(self, iter, base_element, tracers):
        check_conc = self.conc[:,iter]
        # Calculate bgc rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            
            if reac["type"] == "remineralization":  self.remineralization(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)

        if iter % 50 == 0:
            x=1
        x=1
        
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

        tracers[consumed].d_dt -= ec * remineralization
        if "o2" in c:
            if isinstance(parameters["convert_o2"],(int,float)) and not isinstance(parameters["convert_o2"],bool):
                tracers["o2"].d_dt -= remineralization * parameters["convert_o2"]
            elif isinstance(parameters["convert_o2"],str):
                tracers["o2"].d_dt -= remineralization * float(Fraction(parameters["convert_o2"]))
        
        if p[0] == None:    pass
        else:
            tracers[p].d_dt += np.array(ep) * remineralization

    
    