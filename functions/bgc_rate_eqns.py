import numpy as np
import os
import yaml
from functions import rates, seasonal_cycling
from functions.other_functions import concentration_ratio, light_attenuation
from setup.initialize import coordinate_system
np.set_printoptions(precision=20)


# def bgc_rate_eqns(iter, base_element, conc, d_dt, light_attenuation_water, temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind, tracer_map, tracer_type, tracers, sinking):
def bgc_rate_eqns(iter, configuration, base_element, conc, d_dt, light_attenuation_water, temp, sal, dens, z, dz, surface_PAR, wind, tracer_map, tracer_type, tracers, sinking):
    
    # Update concentration ratios
    # for key in tracers:
    #     if tracers[key].type in ["bacteria", "detritus","phytoplankton","zooplankton"]:
    #         # Get index of base element
    #         index = tracers[key].composition.index(base_element)

    #         # Calculate concentration ratios
    #         concentration_ratio(iter, index, tracers[key])

    # if physical["environment"]["forcing"] == "seasonal":
    #     # Seasonal cycling
    #     physical["bgc_phys_vars"]["temperature"] = seasonal_cycling.get_temperature(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_temp"], physical["environment"]["seasonal_cycling"]["summer_temp"])
    #     physical["bgc_phys_vars"]["surface_PAR"] = seasonal_cycling.get_sunlight(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_sun"], physical["environment"]["seasonal_cycling"]["summer_sun"],physical["environment"]["seasonal_cycling"]["latitude"])
    #     physical["bgc_phys_vars"]["mld"] = seasonal_cycling.get_mixed_layer_depth(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_mld"], physical["environment"]["seasonal_cycling"]["summer_mld"])
    #     physical["bgc_phys_vars"]["salinity"] = seasonal_cycling.get_salinity(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_salt"], physical["environment"]["seasonal_cycling"]["summer_salt"])
    #     physical["bgc_phys_vars"]["wind"] = seasonal_cycling.get_wind(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_wind"], physical["environment"]["seasonal_cycling"]["summer_wind"])
    #     physical["bgc_phys_vars"]["z"] = physical["vertical_grid"]["z"]

    # Calculate concentration ratios
    conc_ratio = concentration_ratio(conc, tracer_map)

    # Calculate light attenuation coefficient
    k_PAR = light_attenuation(base_element, light_attenuation_water, conc, tracer_map, tracers)

    # Initialize limitation factor for denitrification to 0.
    bact_limitation_factor = np.zeros_like(z)

    # Calculate bacteria rates (need to do this first to get bact_limitation factor for denitrification)
    for key in tracers:
        if tracers[key].type == "bacteria":
            bact_limitation_factor += tracers[key].bac(base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)
    
    # Calculate other bgc rates
    for key in tracers:
        # if tracers[key].type == "bacteria":
        #     bact_limitation_factor += tracers[key].bac(base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)
        if tracers[key].type == "detritus":
            tracers[key].detritus(base_element, temp, conc, d_dt, tracer_map, tracers)
        elif tracers[key].type == "inorganic": 
            tracers[key].inorg(configuration, bact_limitation_factor, conc, d_dt, tracer_map, z, dz, temp, sal, dens, wind)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(iter, base_element, temp, z, dz, k_PAR, surface_PAR, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers, sinking)
        elif tracers[key].type == "zooplankton": 
            tracers[key].zoo(iter, base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)

    # Convert rates fro 1/d to 1/s
    d_dt /= 86400.

    return d_dt
