import numpy as np

seconds_per_day = 86400
seconds_per_hour = 3600

def light_attenuation(parameters, phyto):

    k_PAR = parameters["light_attenuation_water"] + parameters["light_attenuation_phyto"] * phyto

    return k_PAR

def light_limitation(parameters, irrad, k_PAR, mixed_layer_depth, surface_PAR):
    """
    k_PAR = Light Attenuation Coefficient
    surface_PAR = Photosynthetically Active Radiation (PAR) at watetr surface (z = 0)
    """
    
    # coeff = (parameters["max_photosynthetic_rate"]/seconds_per_hour) / ( k_PAR * mixed_layer_depth )
    # numerator = (parameters["initial_PI_slope"]/seconds_per_hour) * surface_PAR + np.sqrt(((parameters["max_photosynthetic_rate"]/seconds_per_hour)**2) + (((parameters["initial_PI_slope"]/seconds_per_hour) * surface_PAR)**2))
    # denominator = (parameters["initial_PI_slope"]/seconds_per_hour) * irrad + np.sqrt(((parameters["max_photosynthetic_rate"]/seconds_per_hour)**2) + (((parameters["initial_PI_slope"]/seconds_per_hour) * irrad)**2))
    
    coeff = parameters["max_photosynthetic_rate"] / ( k_PAR * mixed_layer_depth )
    numerator = ( parameters["initial_PI_slope"] * surface_PAR ) + np.sqrt( ( parameters["max_photosynthetic_rate"] ** 2 ) + ( ( parameters["initial_PI_slope"] * surface_PAR) ** 2 ) )
    denominator = ( parameters["initial_PI_slope"] * irrad ) + np.sqrt( ( parameters["max_photosynthetic_rate"] ** 2 ) + ( ( parameters["initial_PI_slope"] * irrad ) ** 2 ) )
    
    light_limitation = coeff * np.log(numerator/denominator)

    return light_limitation


def irradiance(surface_PAR, depth, k_PAR):

    irradiance = surface_PAR * np.exp( -k_PAR * depth)

    return irradiance


def max_growth_rate(parameters, temperature):
    """
    Defiition:: Calculates the temperature-dependent maximum phytoplankton grwoth rate, Vm
    """

    Vm = parameters["a"] * ( parameters["b"] ** ( parameters["c"] * temperature ) )

    return Vm


def nutirent_limitation(nutrient, half_saturation_constant):
    """
    Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
    """

    nutrient_limitation_factor = nutrient / (half_saturation_constant + nutrient)
    
    return nutrient_limitation_factor


def temperature_regulation(parameters, base_temp, temperature):

    temp_regulating_factor = np.exp( np.log(parameters["q10_coefficient"]) * (temperature - base_temp) / base_temp )

    return temp_regulating_factor