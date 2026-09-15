import numpy as np
from functions.other_functions import concentration_ratio, light_attenuation
from functions.seasonal_cycling import get_mixed_layer_depth, get_salinity, get_sunlight, get_temperature, get_wind
from pom.calculations import density_profile
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


def reduced_bgc_rate_eqns(time, base_element, conc, tracer_map, tracer_type, tracers, physical):

    # Extract physical variables
    num_layers = physical["water_column"]["num_layers"]                 # number of layers in water column [-]
    column_depth = physical["water_column"]["column_depth"]             # water column depth [m]
    z = physical["vertical_grid"]["z"]                                  # vertical grid [m]
    dz = physical["vertical_grid"]["dz"]                                # vertical spacing [m]
    light_attenuation_water = physical["environment"]["light_attenuation_water"]       # light attenuation coefficient for water

    configuration = physical["simulation"]["configuration"]
    forcing = physical["environment"]["forcing"]
    forcing_data = physical["environment"]["forcing_data"]

    # Define iteration for phytoplankton net primary production storage
    iter = int(np.floor(time/physical["simulation"]["dt"]))
   
    # Create arrays for temperature and salinity if 1D simulation
    if configuration == "1d":
        # Initialize physical variables
        temperature = np.zeros(num_layers-1,dtype=np.float64)
        salinity = np.zeros(num_layers-1,dtype=np.float64)
        mixed_layer_depth = np.zeros(num_layers-1,dtype=np.float64)
        surface_PAR = 0.
        wind = 0.
    
        # Calculate physical variables at current time
        if forcing == "constant":
            temp = forcing_data["temperature"] * np.ones(num_layers-1)
            sal = forcing_data["salinity"] * np.ones(num_layers-1)
            surface_PAR = forcing_data["sunlight"]
            wind = forcing_data["wind"]

        elif forcing == "seasonal":
            # Temperature
            t_win = np.linspace(forcing_data["winter_temp"], 0.8*forcing_data["winter_temp"], num_layers-1)     # (surface_value, bottom_value, steps)
            t_sum = np.linspace(forcing_data["summer_temp"], 0.8*forcing_data["summer_temp"], num_layers-1)     # (surface_value, bottom_value, steps)
            temp = get_temperature(time, t_win, t_sum, forcing_data["temp_excursion"])
            
            # Salinity
            s_win = np.linspace(0.95*forcing_data["winter_salt"], forcing_data["winter_salt"], num_layers-1)    # (surface_value, bottom_value, steps)
            s_sum = np.linspace(0.95*forcing_data["summer_temp"], forcing_data["summer_salt"], num_layers-1)    # (surface_value, bottom_value, steps)
            sal = get_salinity(time, s_win, s_sum)

            # Shortwave irradiance flux
            surface_PAR = get_sunlight(time,forcing_data["winter_sun"], forcing_data["summer_sun"], physical["environment"]["latitude"])

            # Wind speed
            wind = get_wind(time, forcing_data["winter_wind"], forcing_data["summer_wind"])

        dens = density_profile(configuration, num_layers-1, column_depth/2, 0., temp, sal)     # Calculate density in center of cell (column_depth/2)

    elif configuration == "0d":
        # Initialize physical variables
        temperature = np.zeros(num_layers,dtype=np.float64)
        salinity = np.zeros(num_layers,dtype=np.float64)
        mixed_layer_depth = np.zeros(num_layers,dtype=np.float64)
        surface_PAR = 0.
        wind = 0.

        # Calculate physical variables at current time
        if forcing == "constant":
            temperature[0] = forcing_data["temperature"]
            salinity[0] = forcing_data["salinity"]
            mixed_layer_depth[0] = forcing_data["mld"]
            surface_PAR = forcing_data["sunlight"]
            wind[0] = forcing_data["wind"]

        elif forcing == "seasonal":
            temperature[0] = get_temperature(time, forcing_data["winter_temp"], forcing_data["summer_temp"], forcing_data["temp_excursion"])
            salinity[0] = get_salinity(time, forcing_data["winter_salt"], forcing_data["summer_salt"])
            mixed_layer_depth[0] = get_mixed_layer_depth(time,forcing_data["winter_mld"], forcing_data["summer_mld"])
            surface_PAR = get_sunlight(time,forcing_data["winter_sun"], forcing_data["summer_sun"], physical["environment"]["latitude"])
            wind[0] = get_wind(time, forcing_data["winter_wind"], forcing_data["summer_wind"])

        dens = density_profile(configuration, num_layers, column_depth/2, 0., temperature, salinity)     # Calculate density in center of cell (column_depth/2)

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