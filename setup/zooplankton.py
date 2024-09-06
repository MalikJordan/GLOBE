import os
import sys
import numpy as np
from functions.other_functions import concentration_ratio, tracer_elements

class Zooplankton():
    """
    """

    def __init__(self, long_name, composition):
        self.name = long_name
        self.prey = []
        self.grazing_rates = []
        self.composition = []
        self.conc = []
        
        if len(composition) < 1:
            sys.exit("Zooplankton: Element required for " + long_name + ". Check documentation adn edit input file.")
        else:
            for key in composition:
                available_elements = ['c','n','p','fe']
                if key in available_elements:
                    self.composition.append(key)
                    self.conc.append(np.array(composition[key]))
                else:
                    sys.exit("Zooplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
        self.conc = np.array(self.conc)
        self.d_dt = np.zeros_like(self.conc)
        self.conc_ratio = np.zeros_like(self.conc)


    def add_prey(self, prey):
        self.prey.append(prey)


    def egestion(self, parameters, ec, ep, ic, ip, consumed, produced):
        """
        Definition:: Calculates egestion of zooplankton to detrital pool
        Return:: Egestion rate
        """

        egestion = ( 1 - parameters["assimilation_efficiency"] ) * np.sum(self.grazing_rates,axis=0)

        concentration_ratio(ic,consumed)
        concentration_ratio(ip,produced)

        consumed.d_dt -= ec * consumed.conc_ratio * egestion
        produced.d_dt += ep * produced.conc_ratio * egestion


    def excretion(self, parameters, ec, ep, ic, ip, consumed, produced):
        """
        Definition:: Calculates excretion of zooplankton to inorganic nutrient pool.
                     Excretion can be represented either as a constant rate or as a fraction of zooplankton grazing on phytoplankton.
        Return:: Excretion rate
        """

        if parameters["function"] == "constant":
            excretion = parameters["excretion_rate"] * np.array(consumed.conc[ic])
            
        elif parameters["function"] == "grazing":
            excretion = parameters["assimilation_efficiency"] * ( 1 - parameters["ingestion_efficiency"] ) * np.sum(self.grazing_rates,axis=0)

        concentration_ratio(ic,consumed)
        concentration_ratio(ip,produced)

        consumed.d_dt -= ec * consumed.conc_ratio * excretion
        produced.d_dt += ep * produced.conc_ratio * excretion


    def grazing(self, parameters, ec, ep, ic, ip, consumed, produced):
        """
        Definition:: Calculates the zooplankton grazing rate on a particular species using user choice of the Ivlev Equation (1),
                    Holling Type II Response (2), Holling Type III Response (3), or Monod Function (4).
        Return:: grazing_consumed - Grazing rate to be added to "grazing_rates[]" for later summation
                grazing_produced - Porition of grazing rate allocated to zooplaknton (scaled by assimilation and ingestion efficiencies).
        """

        if parameters["function"] == "ivlev":
            function = parameters["max_grazing_rate"] * ( 1 - np.exp( -parameters["ivlev"] * np.array(consumed.conc[ic]) ) )

        elif parameters["function"] == "holling-2":
            function = ( parameters["max_grazing_rate"] * np.array(consumed.conc[ic]) ) / ( parameters["half_sat_grazing"] + np.array(consumed.conc[ic]) )

        elif parameters["function"] == "holling-3":
            function = ( parameters["max_grazing_rate"] * ( np.array(consumed.conc[ic])**2 ) ) / ( (parameters["half_sat_grazing"]**2) + ( np.array(consumed.conc[ic])**2 ) )

        elif parameters["function"] == "monod":
            function = parameters["max_grazing_rate"] * np.array(consumed.conc[ic]) / parameters["half_sat_grazing"]

        grazing_prey = function * np.array(produced.conc[ip])
        grazing_predator = parameters["assimilation_efficiency"] * parameters["ingestion_efficiency"] * function * np.array(produced.conc[ip])

        concentration_ratio(ic,consumed)
        concentration_ratio(ip,produced)

        consumed.d_dt -= ec * consumed.conc_ratio * grazing_prey
        produced.d_dt += ep * produced.conc_ratio * grazing_predator
        produced.grazing_rates.append(grazing_prey)


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


    def respiration():
        pass


class MesoZooplankton(Zooplankton):
    """
    """

    def __init__(self, long_name, composition):
        super().__init__(long_name, composition)
        self.type = "mesozooplankton"

    def mesozoo(self, base_element, reactions, tracers):

        # Zero out grazing rates
        self.grazing_rates = []

        # Calculate bgc rates
        for reac in reactions:
            tc = tracers[reac["consumed"]]  # consumed tracer
            tp = tracers[reac["produced"]]  # produced tracer
            
            # ic/ip = index of base element in consumed/produced tracer
            # ec/ep = array of  elements in consumed/produced tracer affected by current reaction
            ec, ep, ic, ip = tracer_elements(base_element, reac, tc, tp)

            if reac["type"] == "grazing":   super().grazing(reac["parameters"], ec, ep, ic, ip, tc, tp)
            if reac["type"] == "egestion":  super().egestion(reac["parameters"], ec, ep, ic, ip, tc, tp)
            if reac["type"] == "excretion": super().excretion(reac["parameters"], ec, ep, ic, ip, tc, tp)
            if reac["type"] == "mortality": super().mortality(reac["parameters"], ec, ep, ic, ip, tc, tp)


class MicroZooplankton(Zooplankton):
    """
    """

    def __init__(self, long_name, composition):
        super().__init__(long_name, composition)
        self.type = "microzooplankton"
    