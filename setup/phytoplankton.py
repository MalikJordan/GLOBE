import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutirent_limitation, temperature_regulation, tracer_elements

class Phytoplankton():
    """
    
    """

    def __init__(self, long_name, composition, nutrient_limitation):
        self.name = long_name
        self.nutrients = []
        self.nutrient_half_sat = []
        self.nutrient_limitation = nutrient_limitation
        self.nutrient_limitation_factor = 0.
        self.composition = []
        self.conc = []
        self.cff = "living organic"
        self.type = "phytoplankton"

        if len(composition) < 1:
            sys.exit("Phytoplankton: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition:
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    self.composition.append(key)
                    self.conc.append(composition[key])
                else:
                    sys.exit("Phytoplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")

        self.d_dt = np.zeros_like(self.conc)
        self.conc_ratio = np.zeros_like(self.conc)


    def add_nutrient(self, nutrient,half_sat):
        self.nutrients.append(nutrient)
        self.nutrient_half_sat.append(half_sat)

    def calculate_nutrient_limitation(self, tracers):
        """
        Calculates nutrient limitation factor as either a minimum, product, or sum of all nutrients which limit phytoplankton growth
        """
        nutrient_conc = np.zeros(len(self.nutrients))
        for key in tracers:
            if key in self.nutrients:
                i = self.nutrients.index(key)
                nutrient_conc[i] = np.array(tracers[key].conc)
        lim = nutirent_limitation(nutrient_conc, self.nutrient_half_sat)

        if self.nutrient_limitation == "minimum":
            self.nutrient_limitation_factor = np.min(lim, axis=0)
        elif self.nutrient_limitation == "product":
            self.nutrient_limitation_factor = np.prod(lim, axis=0)
        elif self.nutrient_limitation == "sum":
            self.nutrient_limitation_factor == np.sum(lim,  axis=0)


    def mortality(self, parameters, ec, ep, ic, ip, consumed, produced):
        """
        Definition:: This function applies a choice of a linear or quadratic term to account for the non-grazing mortality of planktoninc species.
        Return:: Mortality rate
        """

        if parameters["function"] == "linear":
            mortality = parameters["mortality_rate"] * np.array(consumed.conc[ic])

        elif parameters["function"] == "quadratic":
            mortality = parameters["mortality_rate"] * ( np.array(consumed.conc[ic])**2 )

        concentration_ratio(ic,consumed)
        concentration_ratio(ip,produced)

        consumed.d_dt -= ec * consumed.conc_ratio * mortality
        produced.d_dt += ep * produced.conc_ratio * mortality


    def photosynthesis():
        pass

    def uptake(self, parameters, ep, ip, consumed, produced, base_temp, coordinates, mixed_layer_depth, surface_PAR, temperature):
        
        if parameters["light_limitation"] == "variable":
            k_PAR = light_attenuation(parameters, np.array(produced.conc[ip]))
            irrad = irradiance(surface_PAR, coordinates, k_PAR)
            light_lim = light_limitation(parameters, irrad, k_PAR, mixed_layer_depth, surface_PAR)
        else:
            light_lim = parameters["light_limitation"]

        if parameters["max_growth_rate"] == "variable":
            Vm = max_growth_rate(parameters, temperature)
        else:
            Vm = parameters["max_growth_rate"]
        
        growth = Vm * light_lim * self.nutrient_limitation_factor * np.array(produced.conc[ip])

        concentration_ratio(ip,produced)

        consumed.d_dt -= growth
        produced.d_dt += ep * produced.conc_ratio * growth
        

    def phyto(self, base_element, base_temp, coordinates, mixed_layer_depth, reactions, surface_PAR, temperature, tracers):

        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(tracers)
        
        for reac in reactions:
            tc = tracers[reac["consumed"]]  # consumed tracer
            tp = tracers[reac["produced"]]  # produced tracer
            
            # ic/ip = index of base element in consumed/produced tracer
            # ec/ep = array of  elements in consumed/produced tracer affected by current reaction
            ec, ep, ic, ip = tracer_elements(base_element, reac, tc, tp)

            if reac["type"] == "mortality": self.mortality(reac["parameters"], ec, ep, ic, ip, tc, tp)
            if reac["type"] == "uptake":    self.uptake(reac["parameters"], ep, ip, tc, tp, base_temp, coordinates, mixed_layer_depth, surface_PAR, temperature)
