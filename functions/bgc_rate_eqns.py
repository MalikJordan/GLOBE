import numpy as np
import os
import yaml
from functions import rates, seasonal_cycling
from functions.other_functions import concentration_ratio
from setup.initialize import coordinate_system


def bgc_rate_eqns(iter, base_element, physical, pom1d, tracers):

    # Update concentration ratios
    for key in tracers:
        if tracers[key].type in ["bacteria", "detritus","phytoplankton","zooplankton"]:
            # Get index of base element
            index = tracers[key].composition.index(base_element)

            # Calculate concentration ratios
            concentration_ratio(iter, index, tracers[key])

    if physical["environment"]["forcing"] == "seasonal":
        # Seasonal cycling
        physical["bgc_phys_vars"]["temperature"] = seasonal_cycling.get_temperature(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_temp"], physical["environment"]["seasonal_cycling"]["summer_temp"])
        physical["bgc_phys_vars"]["surface_PAR"] = seasonal_cycling.get_sunlight(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_sun"], physical["environment"]["seasonal_cycling"]["summer_sun"],physical["environment"]["seasonal_cycling"]["latitude"])
        physical["bgc_phys_vars"]["mld"] = seasonal_cycling.get_mixed_layer_depth(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_mld"], physical["environment"]["seasonal_cycling"]["summer_mld"])
        physical["bgc_phys_vars"]["salinity"] = seasonal_cycling.get_salinity(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_salt"], physical["environment"]["seasonal_cycling"]["summer_salt"])
        physical["bgc_phys_vars"]["wind"] = seasonal_cycling.get_wind(physical["simulation"]["time"][iter], physical["environment"]["seasonal_cycling"]["winter_wind"], physical["environment"]["seasonal_cycling"]["summer_wind"])
        physical["bgc_phys_vars"]["z"] = physical["vertical_grid"]["z"]

    # Clear previous rates
    for key in tracers:
        tracers[key].d_dt = np.zeros_like(tracers[key].d_dt)

    # Calculate bgc rates
    for key in tracers:
        if tracers[key].type == "bacteria":
            tracers[key].bac(iter, base_element, physical, tracers)
        elif tracers[key].type == "detritus":
            tracers[key].detritus(iter, base_element, physical, tracers)
        elif tracers[key].type == "inorganic":
            tracers[key].inorg(iter, base_element, physical, tracers)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(iter, base_element, physical, tracers)
        elif tracers[key].type == "zooplankton":
            tracers[key].zoo(iter, base_element, physical, tracers)

    # Convert rates fro 1/d to 1/s
    for key in tracers:
        tracers[key].d_dt /= physical["simulation"]["sec_per_day"]
        
    # # Apply rates to tracer concentrations
    # for key in tracers:
    #     tracers[key].conc[...,iter+1] = tracers[key].conc[...,iter] + physical["simulation"]["dt"] * tracers[key].d_dt / physical["simulation"]["sec_per_day"]  # Convert rates from 1/d to 1/s
    
    # # Set minimum of zero for tracer concentrations
    # for key in tracers:
    #     tracers[key].conc[...,iter+1] = np.maximum(1E-20*np.ones_like(tracers[key].conc[...,iter+1]),tracers[key].conc[...,iter+1])
    
    # # Update concentration ratios
    # for key in tracers:
    #     if tracers[key].type in ["bacteria", "detritus","phytoplankton","zooplankton"]:
    #         # Get index of base element
    #         index = tracers[key].composition.index(base_element)

    #         # Calculate concentration ratios
    #         concentration_ratio(iter+1, index, tracers[key])
