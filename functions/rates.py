import numpy as np
from functions.other_functions import light_attenuation, light_limitation, irradiance, max_growth_rate, nutirent_limitation, temperature_regulation


def egestion(parameters, grazing_rates):

    egestion = ( 1 - parameters["assimilation_efficiency"] ) * np.sum(grazing_rates)

    return egestion


def excretion(parameters, grazing_rates):

    excretion = parameters["assimilation_efficiency"] * ( 1 - parameters["ingestion_efficiency"] ) * np.sum(grazing_rates)

    return excretion


def grazing(parameters, consumed, produced):
    """
    Definition:: Calculates the zooplankton grazing rate on a particular species using user choice of the Ivlev Equation (1),
                 Holling Type II Response (2), or Holling Type III Response (3).
    Return:: Zooplankton grazing rate on tracer

    Maximum Grazing Rate for a particular zooplankton group is the same for all ingested species (i.e. phytoplankotn and/or detritus) 
    but may be different between zooplankton groups.
    """
    if parameters["function"] == "ivlev":
        # Ivlev Equation
        grazing = produced * parameters["max_grazing_rate"] * ( 1 - np.exp( -parameters["ivlev"] * consumed))

    elif parameters["function"] == "holling-2":
        # Holling Type II Response
        # half_sat = parameters["max_grazing_rate"] / parameters["capture_efficiency"]
        # grazing = produced * ( parameters["max_grazing_rate"] * consumed ) / ( half_sat + consumed )
        grazing = produced * ( parameters["max_grazing_rate"] * parameters["capture_efficiency"] * consumed ) / ( parameters["max_grazing_rate"] + ( parameters["capture_efficiency"] * consumed ) )

    elif parameters["function"] == "holling-3":
        # Holling Type III Response
        # half_sat = np.sqrt( parameters["max_grazing_rate"] / parameters["capture_efficiency"] )
        # grazing = produced * ( parameters["max_grazing_rate"] * (consumed**2) ) / ( (half_sat**2) + (consumed**2) )
        grazing = produced * ( parameters["max_grazing_rate"] * parameters["capture_efficiency"] * (consumed**2) ) / ( parameters["max_grazing_rate"] + ( parameters["capture_efficiency"] * (consumed**2) ) )

    return grazing


def mortality(parameters, consumed):
    """
    Definition:: This function applies both linear and non-linear terms to account for the non-grazing mortality of planktoninc species.
    Return:: Natural and total mortality rates

    The linear term may account for metabolic losses or natural mortality.
    The non-linear term is used to represent density-dependent losses (i.e. infections).

    For phytoplankton, both natural and density-dependent losses are allocated to detritus.
    For zooplankton, only natural losses are allocated to detritus. Density dependent losses are considered to be due to predation by carnivores.
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

    Formulation uses q_10 coefficient in place of maximum photosynthetic rate
    """
    
    nut_lim = nutirent_limitation(nutrient, parameters["half_sat_nutrient"])
    
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