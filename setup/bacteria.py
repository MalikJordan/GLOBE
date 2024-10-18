import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Bacteria():
    """
    
    """

    def __init__(self, abbrev, iters, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

         # Nutrient limitation
        self.nutrients = []
        self.nutrient_half_sat = []
        self.nutrient_limitation_type = tracer["parameters"]["nutrient_limitation"]
        self.nutrient_limitation_factor = 0.

        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Bacteria: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p']
                if key in available_elements:
                    self.composition.append(key)
                    conc.append(tracer["composition"][key])
                else:
                    sys.exit("Bacteria: Element '" + key + "' not recognized. Check documentation and edit input file.")
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
            if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
                self.reactions.append(reac)
        
        # Reorder reactions (uptake needs to appear first)
        self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]


    # def __init__(self, abbrev, composition, iters, long_name, parameters, reactions, type):
    #     self.name = long_name
    #     self.type = type

    #      # Nutrient limitation
    #     self.nutrients = []
    #     self.nutrient_half_sat = []
    #     self.nutrient_limitation_type = parameters["nutrient_limitation"]
    #     self.nutrient_limitation_factor = 0.

    #     # Composition and concentration arrays
    #     self.composition = []
    #     conc = []
    #     if len(composition) < 1:
    #         sys.exit("Bacteria: Element required for " + long_name + ". Check documentation adn edit input file.")
    #     else:
    #         for key in composition:
    #             available_elements = ['c','n','p']
    #             if key in available_elements:
    #                 self.composition.append(key)
    #                 conc.append(composition[key])
    #             else:
    #                 sys.exit("Bacteria: Element '" + key + "' not recognized. Check documentation and edit input file.")
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
        
    #     # Reorder reactions (uptake needs to appear first)
    #     self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]


    def add_nutrient(self, nutrient,half_sat):
        self.nutrients.append(nutrient)
        self.nutrient_half_sat.append(half_sat)
        
    def mortality():
        pass

    def uptake():
        pass

    def bac():
        pass
