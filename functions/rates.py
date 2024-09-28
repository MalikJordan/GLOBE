import numpy as np
from functions.other_functions import light_attenuation, light_limitation, max_growth_rate, irradiance, nutrient_limitation, temperature_regulation


def egestion(parameters, grazing_rates):
    """
    Definition:: Calculates egestion of zooplankton to detrital pool
    Return:: Egestion rate
    """

    egestion = ( 1 - parameters["assimilation_efficiency"] ) * np.sum(grazing_rates)

    return egestion


def excretion(parameters, grazing_rates, consumed):
    """
    Definition:: Calculates excretion of zooplankton to inorganic nutrient pool.
                 Excretion can be represented either as a constant rate or as a fraction of zooplankton grazing on phytoplankton.
    Return:: Excretion rate
    """

    if parameters["function"] == "constant":
        excretion = parameters["excretion_rate"] * consumed
        
    elif parameters["function"] == "grazing":
        excretion = parameters["assimilation_efficiency"] * ( 1 - parameters["ingestion_efficiency"] ) * np.sum(grazing_rates)

    return excretion


def grazing(parameters, consumed, produced):
    """
    Definition:: Calculates the zooplankton grazing rate on a particular species using user choice of the Ivlev Equation (1),
                 Holling Type II Response (2), Holling Type III Response (3), or Monod Function (4).
    Return:: grazing_consumed - Grazing rate removed from consumed tracer (full grazing rate)
             grazing_produced - Porition of grazing rate allocated to zooplaknton (scaled by assimilation and ingestion efficiencies).
    """

    if parameters["function"] == "ivlev":
        function = parameters["max_grazing_rate"] * ( 1 - np.exp( -parameters["ivlev"] * consumed ) )

    elif parameters["function"] == "holling-2":
        function = ( parameters["max_grazing_rate"] * consumed ) / ( parameters["half_sat_grazing"] + consumed )

    elif parameters["function"] == "holling-3":
        function = ( parameters["max_grazing_rate"] * (consumed**2) ) / ( (parameters["half_sat_grazing"]**2) + (consumed**2) )

    elif parameters["function"] == "monod":
        function = parameters["max_grazing_rate"] * consumed / parameters["half_sat_grazing"]

    grazing_consumed = function * produced
    grazing_produced = parameters["assimilation_efficiency"] * parameters["ingestion_efficiency"] * function * produced
    
    return grazing_consumed, grazing_produced


def mortality(parameters, consumed):
    """
    Definition:: This function applies a choice of a linear or quadratic term to account for the non-grazing mortality of planktoninc species.
    Return:: Mortality rate
    """

    if parameters["function"] == "linear":
        mortality = parameters["mortality_rate"] * consumed

    elif parameters["function"] == "quadratic":
        mortality = parameters["mortality_rate"] * (consumed**2)

    return mortality


def remineralization(parameters, consumed):
    """
    Definition:: This function calculates the remineralization of nutrients from detrital pools.
    Return:: Detritus remineralization rate
    """

    remineralization = (parameters["remineralization_rate"]) * consumed

    return remineralization


def uptake(parameters, base_temp, coordinates, mixed_layer_depth, nutrient, phyto, surface_PAR, temperature):
    """
    Definition:: Calculates the growth rate of a particular phytoplankton species based on  nutrient uptake
    """
    
    nut_lim = nutrient_limitation(nutrient, parameters["half_sat_nutrient"])
    
    if parameters["light_limitation"] == "variable":
        k_PAR = light_attenuation(parameters, phyto)
        irrad = irradiance(surface_PAR, coordinates, k_PAR)
        light_lim = light_limitation(parameters, irrad, k_PAR, mixed_layer_depth, surface_PAR)
    else:
        light_lim = parameters["light_limitation"]

    if parameters["max_growth_rate"] == "variable":
        Vm = max_growth_rate(parameters, temperature)
    else:
        Vm = parameters["max_growth_rate"]

    # temp_reg = temperature_regulation(parameters, base_temp, temperature)
    # temp_reg = 1

    # growth = Vm * temp_reg * nut_lim * light_lim * phyto
    growth = Vm * nut_lim * light_lim * phyto

    return growth