import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, tracer_elements

class Inorganic():
    """
    
    """

    def __init__(self, long_name, composition):
        self.name = long_name
        self.conc = []
        self.type = "inorganic"

        if len(composition) > 1:
            sys.exit("Inorganic: Only one element accepted per inorganic nutrient. Check documentation adn edit input file.")
        elif len(composition) < 1:
            sys.exit("Inorganic: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition: 
                self.conc.append(composition[key])
            self.d_dt = np.zeros_like(self.conc)

