import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from pom.calculations import temperature_and_salinity_profiles
from functions.bgc_rate_eqns import bgc_rate_eqns
np.set_printoptions(precision=20)

def pom_bgc_1d(iter, configuration, base_element, light_attenuation_water, temperature, salinity, density, inorganic_suspended_matter, shortwave_radiation,
                w_eddy_velocity, w_gen, wind_speed_zonal, wind_speed_meridional, dif_trac,
                dt2, num_layers, vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, column_depth, 
                nrt_o2, nrt_po4, nrt_no3, nrt_nh4, o2b, no3b, ponb_grad, po4b,
                smoth, umolbgc, nbcbgc, ntp, water_specific_heat_times_density, 
                conc_bwd, conc_cur, sinking, tracer_map, tracer_type, tracers):
    
    # Initialize rate of change array
    d_dt = np.zeros_like(conc_cur,dtype=np.float64)

    # Physical variables for bgc rate equations
    temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind = bgc_physical(temperature, salinity, density, inorganic_suspended_matter, vertical_spacing, column_depth, shortwave_radiation, water_specific_heat_times_density, w_eddy_velocity, w_gen, wind_speed_zonal, wind_speed_meridional)

    # Calculate rate of change
    d_dt = bgc_rate_eqns(iter, configuration, base_element, conc_cur, d_dt, light_attenuation_water, temp, sal, dens, z, dz, surface_PAR, wind, tracer_map, tracer_type, tracers, sinking)
    
    if "o2" in tracers: d_o2surf = tracers["o2"].surf_flux
    else:   d_o2surf = np.float64(0.)

    conc_bwd, conc_cur = vertical_diffusivity(iter, conc_bwd, conc_cur, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
                            nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
                            vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, umolbgc, nbcbgc, ntp, shortwave_radiation, dif_trac)

    return conc_bwd, conc_cur


@njit
def bgc_physical(temperature, salinity, density, inorganic_suspended_matter, vertical_spacing, column_depth, swrad, water_specific_heat_times_density, w_eddy_velocity, w_gen, wsu, wsv):

    temp = temperature[:-1]
    sal = salinity[:-1]
    dens = (density[:-1] * 1.E+03) + 1.E+03
    ism = inorganic_suspended_matter[:]
    z = vertical_spacing[:-1] * column_depth
    dz = vertical_spacing[:-1]
    surface_PAR = -swrad * water_specific_heat_times_density
    weddy = w_eddy_velocity[:]
    wgen = w_gen[:]

    rms_wind = np.sqrt(wsu**2 + wsv**2) * 1.E+03
    wind = np.sqrt(rms_wind/(1.25 * 0.0014))

    return temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind


@njit
def vertical_advection(b_cur, b_bwd, b_fwd, sinking_velocity, num_layers, dzr):
    """"
    Description: Handles the sinking of BFM state variablles. Sinking is treated as downward vertical advection
                 computed with upstream finite differences.
    NOTE: Downward velocities are negative
    """
    # sinking velocity input from vdiff_SOS
    b_cur[-1] = b_cur[-2]
    b_bwd[-1] = b_bwd[-2]
    
    b_fwd[0] = dzr[0] * b_cur[0] * sinking_velocity[1]
    for i in range(1,num_layers-1):
        b_fwd[i] = dzr[i] * (b_cur[i] * sinking_velocity[i + 1] - b_cur[i - 1] * sinking_velocity[i])

    return b_fwd


@njit
def vertical_diffusivity(iter, conc_bwd, conc_cur, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
                         nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
                         z, dz, dzz, dzr, umol, nbc, ntp, swrad, kh):
    """
    Description: Calculates the vertical diffusivity of BFM biochemical components and
                 integrates BFM state variables with Source Splitting (SoS) method
    """
    # Reverse the tracer map
    reverse_map = reverse_tracer_map(tracer_map)

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0

    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.  # to m/s

    # Relaxation velocities
    trelax_o2 = nrt_o2 / 86400.
    trelax_po4 = nrt_po4 / 86400.
    trelax_no3 = nrt_no3 / 86400.
    trelax_nh4 = nrt_nh4

    # Loop over bgc state variables
    for i in range(0, len(conc_cur)):   # i = tracer constituent
        # Zeroing of previous tracer
        b_cur = np.zeros(num_layers, dtype=np.float64)
        b_bwd = np.zeros(num_layers, dtype=np.float64)
        b_fwd = np.zeros(num_layers, dtype=np.float64)
        b_surf = 0.
        b_sflx = 0.
        b_bflx = 0.

        # Load BFM state variable
        b_cur[:-1] = conc_cur[i]
        b_bwd[:-1] = conc_bwd[i]

        b_cur[-1] = b_cur[-2]
        b_bwd[-1] = b_bwd[-2]

        # Calculate tracer sinking velocity
        sinking_velocity = W_ON*wgen + Weddy_ON*weddy

        if reverse_map[i] == 'o2':
            b_sflx = -(d_o2surf / 86400.)
            b_bflx = (conc_cur[i,-1] - o2b) * trelax_o2
        elif reverse_map[i] == 'no3':
            b_sflx = 0.
            b_bflx = (conc_cur[i,-1] - no3b) * trelax_no3
        elif reverse_map[i] == 'nh4':
            b_sflx = 0.
            b_bflx = ponb_grad * trelax_nh4
        elif reverse_map[i] == 'po4':
            b_sflx = 0.
            b_bflx = (conc_cur[i,-1] - po4b) * trelax_po4
        elif reverse_map[i] == 'co2':
            b_sflx = 0.
        elif reverse_map[i] == 'sio4':
            b_sflx = 0.
        
        sinking_velocity[:-1] -= sinking[i] / 86400.
        if tracer_type[i] == "phytoplankton":  
            # Final sink value for phytoplankton
            sinking_velocity[-1] = sinking_velocity[-2]

        if tracer_type[i] == "particulate":
            # Final sink value for particulate detritus
            sinking_velocity[-1] = sinking_velocity[-2]

        # Sinking: upstream vertical advection
        b_fwd = vertical_advection(b_cur, b_bwd, b_fwd, sinking_velocity, num_layers, dzr)
        
        # Source splitting (SoS) leapfrog integration
        for j in range(0,num_layers-1):
            b_fwd[j] = b_bwd[j] + dt2*( (b_fwd[j]/column_depth) + d_dt[i,j] ) #+ tracers[key].d_dt[index,i])
        
        # Compute vertical diffusion and terminate integration (implicit leapfrogging)
        b_fwd, b_surf, b_sflx, b_bflx = temperature_and_salinity_profiles('BGC', dt2, num_layers, column_depth, z, dz, dzz, umol, nbc, ntp, swrad, kh, b_fwd, b_surf, b_sflx, b_bflx)
        
        # Clipping (if needed)
        for j in range(0,num_layers-1):
            b_fwd[j] = max(1.E-20,b_fwd[j])
        
        # Mix the time step and restore time sequence
        conc_bwd[i,:] = b_cur[:-1] + 0.5 * smoth * (b_fwd[:-1] + b_bwd[:-1] - 2.*b_cur[:-1])
        conc_cur[i,:] = b_fwd[:-1]

    return conc_bwd, conc_cur


@njit
def get_tracer_from_index(index, tracer_map):
    """
    Definition: Identifies the tracer key associated with index
    :return: tracer key
    """
    for key,value in tracer_map.items():
        if index in value:
            tracer = key

    return tracer


@njit
def reverse_tracer_map(tracer_map):
    """
    Definition: Reverses tracer map for quicker index lookup
    :return: reversed tracer map
    """
    reverse_map = Dict.empty(key_type=types.int64, value_type=types.unicode_type)

    for tracer,value in tracer_map.items():
        for index in value:
            reverse_map[index] = tracer

    return reverse_map
