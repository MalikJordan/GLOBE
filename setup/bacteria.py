import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Bacteria():
    """
    
    """

    def __init__(self, long_name, composition):
        self.name = long_name
        self.nutrients = []
        self.limiting_nutrient = []
        self.composition = []
        self.conc = []
        self.type = "bacteria"

        if len(composition) < 1:
            sys.exit("Bacteria: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition:
                available_elements = ['c','n','p']
                if key in available_elements:
                    self.composition.append(key)
                    self.conc.append(composition[key])
                else:
                    sys.exit("Bacteria: Element '" + key + "' not recognized. Check documentation and edit input file.")

        self.d_dt = np.zeros_like(self.conc)
        self.conc_ratio = np.zeros_like(self.conc)

    def add_nutrient(self, nutrient):
        self.nutrients.append(nutrient)
        self.limiting_nutrient.append(0.)
        
    def mortality():
        pass

    def uptake():
        pass

    def bac():
        pass
