import numpy as np
import sys

def light_attenuation(parameters, phyto):
    """
    Definition:: Calculates light attenuation factor for photosynthesis
    """
    k_PAR = parameters["light_attenuation_water"] + parameters["light_attenuation_phyto"] * phyto

    return k_PAR


def light_limitation(parameters, irrad, k_PAR, mixed_layer_depth, surface_PAR):
    """
    k_PAR = Light Attenuation Coefficient
    surface_PAR = Photosynthetically Active Radiation (PAR) at watetr surface (z = 0)
    """
    coeff = parameters["max_photo_rate"] / ( k_PAR * mixed_layer_depth )
    numerator = ( parameters["initial_PI_slope"] * surface_PAR ) + np.sqrt( ( parameters["max_photo_rate"] ** 2 ) + ( ( parameters["initial_PI_slope"] * surface_PAR) ** 2 ) )
    denominator = ( parameters["initial_PI_slope"] * irrad ) + np.sqrt( ( parameters["max_photo_rate"] ** 2 ) + ( ( parameters["initial_PI_slope"] * irrad ) ** 2 ) )
    
    light_limitation = coeff * np.log(numerator/denominator)

    return light_limitation


def max_growth_rate(parameters, temperature):
    """
    Defiition:: Calculates the temperature-dependent maximum phytoplankton grwoth rate, Vm
    """
    Vm = parameters["a"] * ( parameters["b"] ** ( parameters["c"] * temperature ) )

    return Vm


def irradiance(surface_PAR, depth, k_PAR):
    """
    Definition:: Calculates usable light for photosynthesis
    """
    irradiance = surface_PAR * np.exp( -k_PAR * depth)

    return irradiance


def nutrient_limitation(nutrient, half_sat):
    """
    Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
    """
    nutrient_limitation_factor = nutrient / (half_sat + nutrient)
    
    return nutrient_limitation_factor


def temperature_regulation(base_temp, temperature, q10):
    """
    Definition:: Calculates temperature regulating factor
    """
    temp_regulating_factor = np.exp( np.log(q10) * (temperature - base_temp) / base_temp )

    return temp_regulating_factor


def concentration_ratio(iter, index, tracer):
    """
    Definition:: Calculates concentration ratio of elements in tracer composition to its base element
    """
    for const in range(0,len(tracer.conc[...,iter])):
        tracer.conc_ratio[const] = tracer.conc[const,iter] / (tracer.conc[index,iter] + 1E-20)


def tracer_elements(base_element, reaction, tracers):
    """
    Definition:: Creates dictionary of tracer elements used for a particular reaction

    c/p   = consumed/produced tracer
    ic/ip = index of base element in consumed/produced tracer
    ec/ep = array of  elements in consumed/produced tracer affected by current reaction
    """
    if "consumed" in reaction and reaction["consumed"] != None:    rc = reaction["consumed"]
    else:   rc = None
    if "produced" in reaction and reaction["produced"] != None:    rp = reaction["produced"]
    else:   rp = None

    
    consumed = {}
    ic = {}
    ec = {}

    if rc != None:
        consumed_tracers = list(rc.keys())
        for c in consumed_tracers:
            consumed[c] = tracers[c]
        
        for key in consumed:
            if consumed[key].type == "inorganic":
                ic[key] = 0
                ec[key] = [1.]
            else:
                ic[key] = consumed[key].composition.index(base_element)
                ec[key] = np.zeros(len(consumed[key].composition))
                for element in consumed[key].composition:
                    i = consumed[key].composition.index(element)
                    if element in reaction["consumed"][key]:  ec[key][i] = 1.
    
    produced = {}
    ip = {}
    ep = {}

    if rp != None:
        produced_tracers = list(rp.keys())
        for p in produced_tracers:
            produced[p] = tracers[p]

        for key in produced:
            if produced[key].type == "inorganic":
                ip[key] = 0
                ep[key] = [1.]
            else:
                ip[key] = produced[key].composition.index(base_element)
                ep[key] = np.zeros(len(produced[key].composition))
                for element in produced[key].composition:
                    i = produced[key].composition.index(element)
                    if element in reaction["produced"][key]:  ep[key][i] = 1.
    else:   produced = {None:None}

    c = list(consumed.keys())
    p = list(produced.keys())

    return c, p, ec, ep, ic, ip