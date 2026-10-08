import numpy as np
from functions.other_functions import concentration_ratio, light_attenuation
from functions.seasonal_cycling import get_mixed_layer_depth, get_salinity, get_sunlight, get_temperature, get_wind, calculate_density
np.set_printoptions(precision=20)

def bgc_rate_eqns(iter, configuration, base_element, conc, d_dt, light_attenuation_water, temp, sal, dens, z, dz, surface_PAR, wind, tracer_map, tracer_type, tracers, sinking):

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
        if tracers[key].type == "detritus":
            tracers[key].detritus(base_element, temp, conc, d_dt, tracer_map, tracers)
        elif tracers[key].type == "inorganic": 
            tracers[key].inorg(configuration, bact_limitation_factor, conc, d_dt, tracer_map, z, dz, temp, sal, dens, wind)
        elif tracers[key].type == "phytoplankton":
            tracers[key].phyto(configuration, iter, base_element, temp, z, dz, k_PAR, surface_PAR, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers, sinking)
        elif tracers[key].type == "zooplankton": 
            tracers[key].zoo(iter, base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)

    # Convert rates fro 1/d to 1/s
    d_dt /= 86400.

    return d_dt


# def reduced_bgc_rate_eqns(time, base_element, conc, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names, dic_matrix):
def reduced_bgc_rate_eqns(time, base_element, conc, num_tracers, physical, tracers, tracer_map, tracer_type, indices_to_retain, removed_tracer_names):
    
    # Extract physical variables
    num_boxes = physical["num_boxes"]           # number of boxes in water column [-]
    column_depth = physical["column_depth"]     # water column depth [m]
    z = physical["z"]                           # vertical grid [m]
    dz = physical["dz"]                         # vertical spacing [m]
    light_attenuation_water = physical["light_attenuation_water"]   # light attenuation coefficient for water

    configuration = physical["configuration"]
    forcing = physical["forcing"]
    forcing_data = physical["forcing_data"]

    # Define iteration for phytoplankton net primary production storage, set to 0 because it doesn't have an effect on reduction output
    # iter = int(np.floor(time/physical["dt"]))
    iter = 0

    # Unravel concentration matrix
    conc = conc.reshape((num_tracers,num_boxes))
   
    # Create arrays for temperature and salinity if 1D simulation
    if configuration == "1d":
        # if dic_matrix:
        #     z = physical["dz"][:-1] * column_depth                          # vertical grid [m]
        #     dz = physical["dz"][:-1]                         # vertical spacing [m]
        # else:
        #     z = physical["dz"] * column_depth                          # vertical grid [m]
        #     dz = physical["dz"]                         # vertical spacing [m]

        z = physical["dz"] * column_depth                          # vertical grid [m]
        dz = physical["dz"]                         # vertical spacing [m]
            
        # Initialize physical variables
        temp = np.zeros(num_boxes,dtype=np.float64)
        sal = np.zeros(num_boxes,dtype=np.float64)
        mld = np.zeros(num_boxes,dtype=np.float64)
        surface_PAR = 0.
        wind = 0.
    
        # Calculate physical variables at current time
        if forcing == "constant":
            temp = forcing_data["temperature"] * np.ones(num_boxes)
            sal = forcing_data["salinity"] * np.ones(num_boxes)
            surface_PAR = forcing_data["sunlight"]
            wind = forcing_data["wind"]

        elif forcing == "seasonal":
            # Temperature
            t_win = np.linspace(forcing_data["winter_temp"], 0.8*forcing_data["winter_temp"], num_boxes)     # (surface_value, bottom_value, steps)
            t_sum = np.linspace(forcing_data["summer_temp"], 0.8*forcing_data["summer_temp"], num_boxes)     # (surface_value, bottom_value, steps)
            temp = get_temperature(time, t_win, t_sum, forcing_data["temp_excursion"])
            
            # Salinity
            s_win = np.linspace(0.95*forcing_data["winter_salt"], forcing_data["winter_salt"], num_boxes)    # (surface_value, bottom_value, steps)
            s_sum = np.linspace(0.95*forcing_data["summer_temp"], forcing_data["summer_salt"], num_boxes)    # (surface_value, bottom_value, steps)
            sal = get_salinity(time, s_win, s_sum)

            # Shortwave irradiance flux
            surface_PAR = get_sunlight(time,forcing_data["winter_sun"], forcing_data["summer_sun"], physical["latitude"])

            # Wind speed
            wind = get_wind(time, forcing_data["winter_wind"], forcing_data["summer_wind"])

        dens = calculate_density(temp, sal, z)

        # Clip physical variables

    elif configuration == "0d":
        # Initialize physical variables
        temp = np.zeros(num_boxes,dtype=np.float64)
        sal = np.zeros(num_boxes,dtype=np.float64)
        mld = np.zeros(num_boxes,dtype=np.float64)
        surface_PAR = 0.
        wind = 0.

        # Calculate physical variables at current time
        if forcing == "constant":
            temp[0] = forcing_data["temperature"]
            sal[0] = forcing_data["salinity"]
            mld[0] = forcing_data["mld"]
            surface_PAR = forcing_data["sunlight"]
            wind = forcing_data["wind"]

        elif forcing == "seasonal":
            temp[0] = get_temperature(time, forcing_data["winter_temp"], forcing_data["summer_temp"], forcing_data["temp_excursion"])
            sal[0] = get_salinity(time, forcing_data["winter_salt"], forcing_data["summer_salt"])
            mld[0] = get_mixed_layer_depth(time,forcing_data["winter_mld"], forcing_data["summer_mld"])
            surface_PAR = get_sunlight(time,forcing_data["winter_sun"], forcing_data["summer_sun"], physical["latitude"])
            wind = get_wind(time, forcing_data["winter_wind"], forcing_data["summer_wind"])

        dens = calculate_density(temp, sal, z)

    # Initialize d_dt and sinking arrays
    d_dt = np.zeros_like(conc)
    sinking = np.zeros_like(conc)

    # Calculate concentration ratios
    conc_ratio = concentration_ratio(conc, tracer_map)
    
    # Calculate light attenuation coefficient
    k_PAR = light_attenuation(base_element, light_attenuation_water, conc, tracer_map, tracers)
    
    # Initialize limitation factor for denitrification to 0.
    bact_limitation_factor = np.zeros_like(z)
    
    # Calculate bacteria rates (need to do this first to get bact_limitation factor for denitrification)
    for key in tracers:
        if tracers[key].type == "bacteria" and key not in removed_tracer_names:
            bact_limitation_factor += tracers[key].bac(base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)
    
    # Calculate other bgc rates
    for key in tracers:
        if tracers[key].type == "detritus" and key not in removed_tracer_names:
            tracers[key].detritus(base_element, temp, conc, d_dt, tracer_map, tracers)
        elif tracers[key].type == "inorganic" and key not in removed_tracer_names: 
            tracers[key].inorg(configuration, bact_limitation_factor, conc, d_dt, tracer_map, z, dz, temp, sal, dens, wind)
        elif tracers[key].type == "phytoplankton" and key not in removed_tracer_names:
            tracers[key].phyto(configuration, iter, base_element, temp, z, dz, k_PAR, surface_PAR, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers, sinking)
        elif tracers[key].type == "zooplankton" and key not in removed_tracer_names: 
            tracers[key].zoo(iter, base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)

    for i in range(0,len(d_dt)):
        if i not in indices_to_retain:  d_dt[i,...] = 0.
    # # Calculate bacteria rates (need to do this first to get bact_limitation factor for denitrification)
    # for key in tracers:
    #     if tracers[key].type == "bacteria":
    #         bact_limitation_factor += tracers[key].bac(base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)

    # # Calculate other bgc rates
    # for key in tracers:
    #     if tracers[key].type == "detritus":
    #         tracers[key].detritus(base_element, temp, conc, d_dt, tracer_map, tracers)
    #     elif tracers[key].type == "inorganic": 
    #         tracers[key].inorg(configuration, bact_limitation_factor, conc, d_dt, tracer_map, z, dz, temp, sal, dens, wind)
    #     elif tracers[key].type == "phytoplankton":
    #         tracers[key].phyto(configuration, iter, base_element, temp, z, dz, k_PAR, surface_PAR, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers, sinking)
    #     elif tracers[key].type == "zooplankton": 
    #         tracers[key].zoo(iter, base_element, temp, conc, conc_ratio, d_dt, tracer_map, tracer_type, tracers)
    
    # Convert rates fro 1/d to 1/s
    d_dt /= 86400.

    # Collapse d_dt matrix
    d_dt = d_dt.ravel()
    
    return d_dt