import numpy as np
import sys

def light_attenuation(abbrev, iter, base_element, light_attenuation_water, tracers):
    """
    Definition:: Calculates light attenuation factor for photosynthesis
    Beer's Law attenuation coefficient
    """
    k_PAR = light_attenuation_water * np.ones_like(tracers[abbrev].conc[0,:,iter])

    for key in tracers:
        if tracers[key].type == "detritus":
            base_index = tracers[key].composition.index(base_element)
            k_PAR += tracers[key].light_attenuation * np.array(tracers[key].conc[base_index,:,iter])
        if tracers[key].type == "phytoplankton":
            # Light attenuation coefficient for phytoplankton is calculated using chl if available
            if "chl" in tracers[key].composition:
                chl_index = tracers[key].composition.index("chl")
                k_PAR += tracers[key].light_attenuation * np.array(tracers[key].conc[chl_index,:,iter])
            else:
                base_index = tracers[key].composition.index(base_element)
                k_PAR += tracers[key].light_attenuation * np.array(tracers[key].conc[base_index,:,iter])

    return k_PAR


# def light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm, pl_pc):
def light_limitation(phyto, iter, parameters, dz, irrad, k_PAR, Vm):
    """
    k_PAR = Light Attenuation Coefficient
    surface_PAR = Photosynthetically Active Radiation (PAR) at watetr surface (z = 0)
    """
    # -------------------------------------------------------------------------------------------------
    # Monod
    # -------------------------------------------------------------------------------------------------
    if parameters["light_limitation"] == "monod":
        # Heinle & Slawig (2013)
        light_limitation = irrad / (parameters["half_sat_light"] + irrad + 1E-20)

    # -------------------------------------------------------------------------------------------------
    # Geider et al. (1997) / Jassby and Platt (1976)
    # -------------------------------------------------------------------------------------------------
    elif parameters["light_limitation"] in ["geider","platt"]:
        
        # Calculate irradiance at depth
        if parameters["light_location"] == "top":
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * 86400                              # *86400 to convert from [1/s] to [1/d]
        elif parameters["light_location"] == "middle":  # Lazzari et al. (2012)
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * np.exp( -k_PAR * dz/2) * 86400     # *86400 to convert from [1/s] to [1/d]
        elif parameters["light_location"] == "integrated":  # Vichi et al. (2007)
            r = irrad / (k_PAR * dz) * (1. - np.exp(-k_PAR*dz))            
            irrad_at_depth = np.maximum(1E-20*np.ones_like(r), r*86400)                                        # *86400 to convert from [1/s] to [1/d]

        # Calculate Chl:C ratio using either ...   
        if {"c","chl"}.issubset(phyto.composition): # Carbon and Chlorophyll concentrations (if preselt)
        # if "c" in self.composition and "chl"  in self.composition:
            carbon_index = phyto.composition.index("c")
            pc = phyto.conc[carbon_index,:,iter]
            chl_index = phyto.composition.index("chl")
            pl = phyto.conc[chl_index,:,iter]
            pl_pc = pl / pc     # Chl:C ratio (used in light limitation)
        else: # Geider et al. (1997) Dynamic Model
            if "theta_min" not in parameters: parameters["theta_min"] = 0.
            pl_pc = (parameters["theta_max"] - parameters["theta_min"]) / ( 1. + ( ( parameters["theta_max"] * parameters["initial_PI_slope"] * irrad_at_depth ) / 
                                                                                ( 2 * Vm * phyto.temp_regulation_factor * phyto.nutrient_limitation_factor + 1.E-20 ) ) ) \
                        + parameters["theta_min"]

        # Calculate exponent for light limitation
        exp = pl_pc * ( parameters["initial_PI_slope"] / Vm ) * irrad_at_depth     # Stays like this for Platt
        if parameters["light_limitation"] == "geider":  # scale by temperature and nutrient limitation factors for Geider
            exp = exp / ( phyto.temp_regulation_factor * phyto.nutrient_limitation_factor )

        light_limitation = 1. - np.exp(-exp)

    # -------------------------------------------------------------------------------------------------
    # Smith (1936)
    # -------------------------------------------------------------------------------------------------
    elif parameters["light_limitation"] == "smith":
        # Evans & Parslow (1985) formulation
        num = Vm * parameters["initial_PI_slope"] * irrad
        den = np.sqrt((Vm**2) + ((parameters["initial_PI_slope"]*irrad)**2))
        
        light_limitation = num/(den + 1E-20)

    return exp, irrad_at_depth, light_limitation


def max_growth_rate(parameters, temperature):
    """
    Defiition:: Calculates the temperature-dependent maximum phytoplankton growth rate, Vm
    """
    if parameters["type"] == "base_b":
        Vm = parameters["a"] * ( parameters["b"] ** ( parameters["c"] * temperature ) )
    elif parameters["type"] == "standard":
        Vm = parameters["base_growth_rate"] * np.exp(parameters["eppley_coeff"] * temperature)

    return Vm


def irradiance(eps_PAR, surface_PAR, depth, k_PAR):
    """
    Definition:: Calculates usable light for photosynthesis
    eps_PAR = fraction of photosynthetically available radiation
    0.217 = conversion from Einstein to Watts
    """

    # irradiance = surface_PAR * eps_PAR / 0.217

    irradiance = np.zeros(len(depth))
    irradiance[0] = surface_PAR * eps_PAR / 0.217
    if len(depth) > 1:
        for i in range(1,len(depth)):
            irradiance[i] = irradiance[i-1] * np.exp(-1. * k_PAR[i-1] * depth[i-1])

    return irradiance


def nutrient_limitation(nutrient, half_sat):
    """
    Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
    """
    nutrient_limitation_factor = nutrient / (half_sat + nutrient + 1.E-20)
    
    return nutrient_limitation_factor


def monod(nutrient, half_sat, exponent):

    limitation_factor = np.power(nutrient, exponent) / ( np.power(nutrient, exponent) + np.power(half_sat, exponent) + 1.E-20)

    return limitation_factor

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


# def temperature_dependence(base_temp, temperature, q10, tracer):
def temperature_dependence(temperature, tracer):
    
    """
    Definition:: Calculates temperature regulating factor
    """
    if tracer.temperature_regulation["function"] == "arrhenius":
        # Convert temperature from Celsius to Kelvin (+273.15)
        # Universal gas constant (R = 8.314 J mol-1 K-1)
        temp_regulating_factor = tracer.temperature_regulation["coefficient"] * np.exp(-tracer.temperature_regulation["activation_energy"] / (8.314 * (temperature+273.15)) )
    
    elif tracer.temperature_regulation["function"] == "eppley":
        temp_regulating_factor = np.exp(tracer.temperature_regulation["coefficient"] * temperature)
    
    elif tracer.temperature_regulation["function"] == "q10":
        # temp_regulating_factor = q10**((temperature-base_temp)/base_temp)
        temp_regulating_factor = tracer.temperature_regulation["coefficient"]**((temperature - tracer.temperature_regulation["base_temp"]) / tracer.temperature_regulation["base_temp"])

    return temp_regulating_factor


def concentration_ratio(iter, index, tracer):
    """
    Definition:: Calculates concentration ratio of elements in tracer composition to its base element
    """
    # for const in range(0,len(tracer.conc[...,iter])):
    #     tracer.conc_ratio[const] = tracer.conc[const,iter] / (tracer.conc[index,iter] + 1E-20)

    for const in range(0,len(tracer.conc)):
        tracer.conc_ratio[const] = tracer.conc[const,:,iter] / (tracer.conc[index,:,iter] + 1E-20)


        # Detritus may be initialized to zero, fix concentration ratio
        # if iter == 0 and tracer.type == 'detritus':
        #     if tracer.conc[const,:,iter].all() == np.zeros(len(tracer.conc[const,:,iter])):
        #         tracer.conc_ratio[const,:] = 1.   # Initialize first concentration ratio for detritus to 1
    
    # Concentration ratio of base element is alway 1
    tracer.conc_ratio[index,:] = 1.


def concentration_ratio_solveivp(concentration, index, indices, tracer):
    """
    Definition:: Calculates concentration ratio of elements in tracer composition to its base element
    """
    for const in range(0,len(indices)):
        tracer.conc_ratio[const] = np.maximum(1E-20, concentration[indices[const]] / (concentration[indices[index]] + 1E-20))
    
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


# def on_off(parameter):

#     x = len(parameter)
#     if 


def switch(parameter):

    x = len(parameter)
    if x > 1:
        x = np.shape(parameter)[0]
        switch = np.zeros(x)
        for i in range(0,x):
            if parameter[i] > 0.0:
                switch[i] = 1.0
    else:
        if parameter > 0.:  switch = 1.
        else:   switch = 0.

    return switch