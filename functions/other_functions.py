import numpy as np
import sys

def light_attenuation(parameters, phyto):
    """
    Definition:: Calculates light attenuation factor for photosynthesis
    Beer's Law attenuation coefficient
    """
    k_PAR = parameters["light_attenuation_water"] + parameters["light_attenuation_phyto"] * phyto

    return k_PAR


# def light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm, pl_pc):
def light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm):
    """
    k_PAR = Light Attenuation Coefficient
    surface_PAR = Photosynthetically Active Radiation (PAR) at watetr surface (z = 0)
    """
    if parameters["light_limitation"] == "monod":
        # Heinle & Slawig (2013)
        light_limitation = irrad / (parameters["half_sat_light"] + irrad + 1E-20)

    elif parameters["light_limitation"] == "platt":  # Jassby and Platt (1976)
        # *86400 to convert from [1/s] to [1/d]
        if parameters["light_location"] == "top":
            r = np.maximum(1E-20*np.ones_like(irrad), irrad) * 86400
        elif parameters["light_location"] == "middle":
            r = np.maximum(1E-20*np.ones_like(irrad), irrad) * np.exp( -k_PAR * dz/2) * 86400
        elif parameters["light_location"] == "integrated":
            r = irrad / (k_PAR * dz) * (1. - np.exp(-k_PAR*dz))
        irr = np.maximum(1E-20*np.ones_like(r), r*86400)    
        exp = pl_pc * ( parameters["max_light_utilization"] / Vm ) * irr

        light_limitation = 1. - np.exp(-exp)

    elif parameters["light_limitation"] == "smith": # Smith (1936)
        # Evans & Parslow (1985) formulation
        num = Vm * parameters["initial_PI_slope"] * irrad
        den = np.sqrt((Vm**2) + ((parameters["initial_PI_slope"]*irrad)**2))
        
        light_limitation = num/(den + 1E-20)

        # Anderson et al. (2015) formulation
        # coeff = Vm / ( k_PAR * mixed_layer_depth )
        # numerator = ( parameters["initial_PI_slope"] * surface_PAR ) + np.sqrt( ( Vm ** 2 ) + ( ( parameters["initial_PI_slope"] * surface_PAR) ** 2 ) )
        # denominator = ( parameters["initial_PI_slope"] * irrad ) + np.sqrt( ( Vm ** 2 ) + ( ( parameters["initial_PI_slope"] * irrad ) ** 2 ) )
        
        # light_limitation = coeff * np.log(numerator/denominator)

    return light_limitation


def max_growth_rate(parameters, temperature):
    """
    Defiition:: Calculates the temperature-dependent maximum phytoplankton grwoth rate, Vm
    """
    Vm = parameters["a"] * ( parameters["b"] ** ( parameters["c"] * temperature ) )

    return Vm


def irradiance(eps_PAR, surface_PAR, depth, k_PAR):
    """
    Definition:: Calculates usable light for photosynthesis
    eps_PAR = fraction of photosynthetically available radiation
    0.217 = conversion from Einstein to Watts
    """
    # irradiance = ( surface_PAR * eps_PAR / 0.217) * np.exp( -k_PAR * depth)
    irradiance = surface_PAR * eps_PAR / 0.217

    return irradiance


def nutrient_limitation(nutrient, half_sat):
    """
    Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
    """
    nutrient_limitation_factor = nutrient / (half_sat + nutrient)
    
    return nutrient_limitation_factor

# def nutrient_limitation(self, tracers):
#     """
#     Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
#     """
#     # Locate index of nutrient element in composition
#     for nut in self.nutrient_limitation[]
#     fN = np.zeros(len(self.nutrient_limitation["nutrients"]),dtype=np.ndarray)
#     if self.nutrient_limitation["type"] == "external":
#         fN = np.minimum(np.ones_like())
    
#     return fN


def temperature_regulation(base_temp, temperature, q10):
    """
    Definition:: Calculates temperature regulating factor
    """
    # temp_regulating_factor = np.exp( np.log(q10) * (temperature - base_temp) / base_temp )
    temp_regulating_factor = q10**((temperature-base_temp)/base_temp)

    return temp_regulating_factor


def concentration_ratio(iter, index, tracer):
    """
    Definition:: Calculates concentration ratio of elements in tracer composition to its base element
    """
    for const in range(0,len(tracer.conc[...,iter])):
        tracer.conc_ratio[const] = tracer.conc[const,iter] / (tracer.conc[index,iter] + 1E-20)
    
    # Concentration ratio of base element is alway 1
    tracer.conc_ratio[index] = 1.


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