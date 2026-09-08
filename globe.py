import os
import time
# import sys
import numpy as np
import yaml
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from setup.initialize import import_bgc_model, import_physical_model
from functions.seasonal_cycling import get_mixed_layer_depth, get_salinity, get_sunlight, get_temperature, get_wind
from functions.bgc_rate_eqns import bgc_rate_eqns
from functions.calculate_averages import average
from pom.calculations import density_profile, kinetic_energy_profile, temperature_and_salinity_profiles, zonal_velocity_profile, meridional_velocity_profile
from pom.forcing import forcing_manager
from pom.initialize import initialize_pom
from pom.coupling import pom_bgc_1d
np.set_printoptions(precision=20)

# def create_function_inputs(tracers):
#     """
#     Definition: Takes tracer dictionary and creates lists, arrays, or typed.Dicts for numba calculations

#     :return: concentration (array), sinking velocities (array), tracer map (typed.Dict), tracer types (list)
#     """

#     # Create list of concentrations
#     concentration = []

#     # Create typed.Dict of tracer indices in concentration
#     tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
    
#     # Create list of trcaer types
#     tracer_type = []   # used in vertical diffusivity calculations

#     # Create list of sinking velocities for each tracer
#     sinking = []

#     index = 0   # counting number for tracer indices
#     for trac in tracers:
#         num_constituents = len(tracers[trac].composition)   # number of constituents in tracer

#         lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
#         for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
#         tracer_map[trac] = lst  # identify tracer constituents with their own index

#         for i in range(num_constituents):
#             # add concentration to matrix
#             concentration.append(tracers[trac].conc[i,...])    # add concentration to matrix

#             # add tracer type to list
#             if tracers[trac].type == "detritus":    tracer_type.append(tracers[trac].form)     # need to distinguish particulate/dissolved form
#             else:   tracer_type.append(tracers[trac].type)     # just the type

#             # add sinking velocity to list
#             if hasattr(tracers[trac],"sinking_velocity"):   sinking.append(tracers[trac].sinking_velocity)
#             else:   sinking.append(np.zeros(tracers[trac].conc.shape[1]))

#             # add tracer type to list
#             index += 1  # update index

#     concentration = np.array(concentration,dtype=np.float64)    # convert concentration from list to array
#     sinking = np.array(sinking,dtype=np.float64)    # convert sinking from list to array

#     return concentration, sinking, tracer_map, tracer_type

def create_function_inputs(iters, tracers):
    """
    Definition: Takes tracer dictionary and creates lists, arrays, or typed.Dicts for numba calculations

    :return: concentration (array), sinking velocities (array), tracer map (typed.Dict), tracer types (list)
    """

    # Create list of concentrations
    initial_concentration = []

    # Create typed.Dict of tracer indices in concentration
    tracer_map = Dict.empty(key_type=types.unicode_type, value_type=types.ListType(types.int64))
    
    # Create list of trcaer types
    tracer_type = []   # used in vertical diffusivity calculations

    # Create list of sinking velocities for each tracer
    sinking = []

    index = 0   # counting number for tracer indices
    for trac in tracers:
        num_constituents = len(tracers[trac].composition)   # number of constituents in tracer

        lst = List.empty_list(types.int64)  # empty typed.List to store elements for tracer constituents
        for i in range(index,index+num_constituents):  lst.append(np.int64(i))  # fill list
        tracer_map[trac] = lst  # identify tracer constituents with their own index

        for i in range(num_constituents):
            # add concentration to matrix
            initial_concentration.append(tracers[trac].initial_conc[i,...])    # add concentration to matrix

            # add tracer type to list
            if tracers[trac].type == "detritus":    tracer_type.append(tracers[trac].form)     # need to distinguish particulate/dissolved form
            else:   tracer_type.append(tracers[trac].type)     # just the type

            # add sinking velocity to list
            if hasattr(tracers[trac],"sinking_velocity"):   sinking.append(tracers[trac].sinking_velocity)
            else:   sinking.append(np.zeros(tracers[trac].initial_conc.shape[1]))

            # add tracer type to list
            index += 1  # update index

    initial_concentration = np.array(initial_concentration,dtype=np.float64)    # convert concentration from list to array
    concentration = np.zeros((initial_concentration.shape[0],initial_concentration.shape[1],iters),dtype=np.float64)
    concentration[:,:,0] = initial_concentration.copy()
    sinking = np.array(sinking,dtype=np.float64)    # convert sinking from list to array

    return concentration, sinking, tracer_map, tracer_type


start = time.perf_counter()
# from pom.check_phys import dens, u, ub, v, vb, t, tb, s, sb, q2, q2b, q2l, q2lb, km, kh, kq
# ----------------------------------------------------------------------------------------------------
# Import and initialize model
# ----------------------------------------------------------------------------------------------------
# check_file = False
# first_check = True
# while not check_file:
#     if first_check:
#         file = input("\nEnter the name of the yaml input file. Press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     else:
#         file = input("\nInput file not found. Re-enter the name of the yaml input file or press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     if file == 'Q' or file == 'q':
#         sys.exit()
#     check_file = os.path.exists(file)
#     first_check = False

# Import physical model
# file = 'physical_bfm17_1d.yaml'
file = 'physical_bfm56.yaml'
file_path = os.getcwd() + '/' + file
physical = import_physical_model(file_path)

# file = 'bfm17_1d.yaml'
file = 'bfm56.yaml'
file_path = os.getcwd() + '/' + file
base_element, reactions, tracers = import_bgc_model(file_path, physical)

concentration, sinking, tracer_map, tracer_type = create_function_inputs(physical["simulation"]["iters"],tracers)

# ----------------------------------------------------------------------------------------------------
# Extract commonly used variables to avoid repetitive dictionary unpacking (unchanged through simulation)
# ----------------------------------------------------------------------------------------------------
iters = physical["simulation"]["iters"]                             # iterations needed
dt = physical["simulation"]["dt"]                                   # time step [s]
dt2  = physical["simulation"]["dt2"]                                # twice the time step [s]
num_layers = physical["water_column"]["num_layers"]                 # number of layers in water column [-]
column_depth = physical["water_column"]["column_depth"]             # water column depth [m]
z = physical["vertical_grid"]["z"]                                  # vertical grid [m]
zz = physical["vertical_grid"]["zz"]                                # staggered vertical grid [m]
dz = physical["vertical_grid"]["dz"]                                # vertical spacing [m]
dzz = physical["vertical_grid"]["dzz"]                              # staggered vertical spacing [m]
dzr = physical["vertical_grid"]["dzr"]                              # reciprocal of vertical spacing [m^-1]
upper_depth = physical["water_column"]["upper_depth"]               # depth before logarithmic spacing
lambda_w = physical["environment"]["light_attenuation_water"]       # light attenuation coefficient for water

configuration = physical["simulation"]["configuration"]

# ----------------------------------------------------------------------------------------------------
# Initialize POM1D (if necessary)
# ----------------------------------------------------------------------------------------------------
# if physical["environment"]["forcing"] == "pom1d":
#     with open(os.getcwd() + '/pom1d.yaml', 'r') as f:
#         pom1d = yaml.full_load(f)
#     physical, forcing = initialize_pom(pom1d, physical)
#     physical = density_profile(physical)
#     # pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * np.pi / 360.)
#     pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * (3.14159265359) / 360.)

if physical["environment"]["forcing"] == "pom1d":
    with open(os.getcwd() + '/pom1d.yaml', 'r') as f:
        pom1d = yaml.full_load(f)

    dif_trac, dif_mom, dif_ke, ke_cur, ke_bwd, ke_fwd, kel_cur, kel_bwd, kel_fwd, \
    u_cur, u_bwd, u_fwd, v_cur, v_bwd, v_fwd, \
    temp_cur, temp_bwd, temp_fwd, temp_adv, temp_int, temp_surf, temp_sflx, temp_bflx, \
    sal_cur, sal_bwd, sal_fwd, sal_adv, sal_int, sal_surf, sal_sflx, sal_bflx, \
    no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb_grad, \
    wsu, wsv, bsu, bsv, ism, swrad, wgen, weddy, mld, \
    counter_ids, counter_params, forcing_ids, forcing_month1, forcing_month2        = initialize_pom(num_layers, pom1d["input_files"]["temperature_IC"], pom1d["input_files"]["salinity_IC"])
    
    # density = density_profile(num_layers, column_depth, physical["vertical_grid"]["dzz"], temp_cur, sal_cur)
    # density = density_profile(num_layers, column_depth, physical["vertical_grid"]["dzz"], temp_bwd, sal_bwd)
    density = density_profile(configuration, num_layers, column_depth, physical["vertical_grid"]["dzz"], temp_bwd, sal_bwd)
    # pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * np.pi / 360.)
    pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * (3.14159265359) / 360.)

    # ----------------------------------------------------------------------------------------------------
    # Extract commonly used variables to avoid repetitive dictionary unpacking (unchanged through simulation)
    # ----------------------------------------------------------------------------------------------------
    idiagn = pom1d["general"]["idiagn"]                             # flag for prognostic/diagnostic mode
    coriolis = pom1d["general"]["coriolis"]                         # coriiolis parameter [s^-1]
    umol = pom1d["background_diffusion"]["umol"]                    # background diffusion coefficient (general)
    umolt = pom1d["background_diffusion"]["umolt"]                  # background diffusion coefficient (temperature)
    umols = pom1d["background_diffusion"]["umols"]                  # background diffusion coefficient (salinity)
    umolbgc = pom1d["background_diffusion"]["umolbgc"]              # background diffusion coefficient (bgc)    
    trt = pom1d["relaxation_times"]["trt"]                          # relaxation time (temperature)
    srt = pom1d["relaxation_times"]["srt"]                          # relaxation time (salinity)
    ssrt = pom1d["relaxation_times"]["ssrt"]                        # relaxation time (surface salinity flux)
    smoth = pom1d["general"]["smoth"]                               # Asselin filter temporal smoother
    nbct = pom1d["flags"]["nbct"]                                   # flag for temperature boundary condition
    nbcs = pom1d["flags"]["nbcs"]                                   # flag for salinity boundary condition
    nbcbgc = pom1d["flags"]["nbcbgc"]                               # flag for bgc boundary condition
    ntp = pom1d["flags"]["ntp"]                                     # flag for jerlov water type
    rcp = pom1d["general"]["water_specific_heat_times_density"]     # specific heat times rho0
    nrt_o2 = pom1d["relaxation_velocities"]["nrt_o2"]               # relaxation velocity for o2
    nrt_po4 = pom1d["relaxation_velocities"]["nrt_po4"]             # relaxation velocity for po4
    nrt_no3 = pom1d["relaxation_velocities"]["nrt_no3"]             # relaxation velocity for no3
    nrt_nh4 = pom1d["relaxation_velocities"]["nrt_nh4"]             # relaxation velocity for nh4

    # Input file strings
    wind_inp = pom1d["input_files"]["wind_stress"]
    rad_inp =  pom1d["input_files"]["shortwave_solar_radiation"]
    heat_inp = pom1d["input_files"]["heat_flux_loss"]
    ism_inp = pom1d["input_files"]["inorganic_suspended_matter"]
    surf_sal_inp = pom1d["input_files"]["surface_salinity"]
    sal_inp = pom1d["input_files"]["salinity"]
    sal_IC_inp = pom1d["input_files"]["salinity_IC"]
    temp_inp = pom1d["input_files"]["temperature"]
    temp_IC_inp = pom1d["input_files"]["temperature_IC"]
    w_vel_inp = pom1d["input_files"]["w_velocity"]
    weddy1_inp = pom1d["input_files"]["eddy_w_velocity_1"]
    weddy2_inp = pom1d["input_files"]["eddy_w_velocity_2"]
    surf_nut_inp = pom1d["input_files"]["surface_nutrients"]
    bot_nut_inp = pom1d["input_files"]["bottom_nutrients"]
    input_files = [wind_inp, rad_inp, heat_inp, ism_inp, surf_sal_inp, sal_inp, sal_IC_inp, temp_inp, temp_IC_inp, w_vel_inp,
                   weddy1_inp, weddy2_inp, surf_nut_inp, bot_nut_inp]
    

# elif physical["environment"]["forcing"] == "seasonal":
#     seasonal = physical["environment"]["seasonal_cycling"]
#     time_array = physical["simulation"]["time"]
    
#     for iter in range(0,iters-1):
#         temperature = get_temperature(time_array[iter], seasonal["winter_temp"], seasonal["summer_temp"], seasonal["temp_excursion"])
#         salinity = get_salinity(time_array[iter], seasonal["winter_salt"], seasonal["summer_salt"])
#         mixed_layer_depth = get_mixed_layer_depth()
#         surfacer_PAR = get_sunlight()
#         wind = get_wind(time_array[iter], seasonal["winter_wind"], seasonal["summer_wind"])

#         density = density_profile(num_layers, column_depth, physical["vertical_grid"]["dzz"], temperature, salinity)

# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
for iter in range(0,iters-1):

    # Turbulence closure
    ke_fwd[:] = ke_bwd[:]
    kel_fwd[:] = kel_bwd[:]

    ke_fwd, kel_fwd, physical["vertical_grid"]["l"], dif_mom, dif_trac, dif_ke \
        = kinetic_energy_profile(dt2, num_layers, column_depth, z, dz, dzz, physical["vertical_grid"]["l"], umol, dif_mom, dif_trac, dif_ke, 
                                 ke_cur, ke_fwd, ke_bwd, kel_cur, kel_fwd, kel_bwd, density, u_cur, bsu, wsu, v_cur, bsv, wsv)

    # Define forcings
    forcing_month1, forcing_month2, temp_fwd, temp_int, temp_surf, temp_sflx, sal_fwd, sal_int, sal_surf, \
        ism, swrad, wsu, wsv, wgen, weddy, no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb_grad \
            = forcing_manager(iter, dt, num_layers, pom1d, counter_params, forcing_month1, forcing_month2)

    # Temperature and salinity computation
    if idiagn == 0: # Prognostic mode
        # Temperature and salinity fully computed by model
        temp_surf = temp_fwd[0]
        sal_surf = sal_fwd[0]
        
        if trt != 0:
            for j in range(0, num_layers):
                if (-dzz[j] * column_depth) >= upper_depth:
                    temp_adv[j] = (temp_int[j] - temp_cur[j]) / (trt * 86400)

        if srt != 0:
            for j in range(0, num_layers):
                if (-dzz[j] * column_depth) >= upper_depth:
                    sal_adv[j] = (sal_int[j] - sal_cur[j]) / (srt * 86400)
        
        # Calculate surface salinity flux
        sal_sflx = -(sal_surf - sal_cur[0]) * ssrt / 86400

        # Calculate temperature
        temp_fwd[:] = temp_bwd[:] + (temp_adv[:] * dt2)
        temp_fwd, temp_surf, temp_sflx, temp_bflx = temperature_and_salinity_profiles('Temperature', dt2, num_layers, column_depth, z, dz, dzz, 
                                                                                      umolt, nbct, ntp, swrad, dif_trac, temp_fwd, temp_surf, temp_sflx, temp_bflx)
        
        # Calculate salinity
        sal_fwd[:] = sal_bwd[:] + (sal_adv[:] * dt2)
        sal_fwd, sal_surf, sal_sflx, sal_bflx = temperature_and_salinity_profiles('Salinity', dt2, num_layers, column_depth, z, dz, dzz, 
                                                                                  umols, nbcs, ntp, swrad, dif_trac, sal_fwd, sal_surf, sal_sflx, sal_bflx)

        # Mix the timestep (Asselin filter)
        temp_cur[:] = temp_cur[:] + 0.5 * smoth * (temp_fwd[:] + temp_bwd[:] - 2. * temp_cur[:])
        sal_cur[:] = sal_cur[:] + 0.5 * smoth * (sal_fwd[:] + sal_bwd[:] - 2. * sal_cur[:])

    # Velocity computation
    u_fwd[:] = u_bwd[:] + dt2 * coriolis * v_cur[:]
    u_fwd, bsu = zonal_velocity_profile(dt2, num_layers, column_depth, dz, dzz, umol, dif_mom, u_fwd, bsu, wsu)

    v_fwd[:] = v_bwd[:] - dt2 * coriolis * u_cur[:]
    v_fwd, bsv = meridional_velocity_profile(dt2, num_layers, column_depth, dz, dzz, umol, dif_mom, v_fwd, bsv, wsv)

    # Mix the timestep (Asselin filter)
    ke_cur[:] = ke_cur[:] + 0.5 * smoth * (ke_fwd[:] + ke_bwd[:] - 2. * ke_cur[:])
    kel_cur[:] = kel_cur[:] + 0.5 * smoth * (kel_fwd[:] + kel_bwd[:] - 2. * kel_cur[:])

    u_cur[:] = u_cur[:] + 0.5 * smoth * (u_fwd[:] + u_bwd[:] - 2. * u_cur[:])
    v_cur[:] = v_cur[:] + 0.5 * smoth * (v_fwd[:] + v_bwd[:] - 2. * v_cur[:])

    # Restore the time sequence
    ke_bwd[:] = ke_cur[:]
    ke_cur[:] = ke_fwd[:]
    kel_bwd[:] = kel_cur[:]
    kel_cur[:] = kel_fwd[:]

    u_bwd[:] = u_cur[:]
    u_cur[:] = u_fwd[:]
    v_bwd[:] = v_cur[:]
    v_cur[:] = v_fwd[:]

    temp_bwd[:] = temp_cur[:]
    temp_cur[:] = temp_fwd[:]
    sal_bwd[:] = sal_cur[:]
    sal_cur[:] = sal_fwd[:]

    # Update density
    # density = density_profile(num_layers, column_depth, dzz, temp_cur, sal_cur)
    # density = density_profile(num_layers, column_depth, dzz, temp_bwd, sal_bwd)
    density = density_profile(configuration, num_layers, column_depth, dzz, temp_cur, sal_cur)
    # density = density_profile(configuration, num_layers, column_depth, dzz, temp_bwd, sal_bwd)

    # if iter < 10:
    #     filename = f"diffusion_iter{iter:01d}.npz"
    #     load_diffusion = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename, allow_pickle=True)

    #     delta_mom = dif_mom - load_diffusion["mom"]
    #     delta_trac = dif_trac - load_diffusion["trac"]
    #     delta_ke = dif_ke - load_diffusion["ke"]

    #     filename2 = f"velocity_iter{iter:01d}.npz"
    #     load_velocity = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename2, allow_pickle=True)

    #     delta_u_cur = u_cur - load_velocity["u_cur"]
    #     delta_u_bwd = u_bwd - load_velocity["u_bwd"]
    #     delta_u_fwd = u_fwd - load_velocity["u_fwd"]
    #     delta_v_cur = v_cur - load_velocity["v_cur"]
    #     delta_v_bwd = v_bwd - load_velocity["v_bwd"]
    #     delta_v_fwd = v_fwd - load_velocity["v_fwd"]

    #     x=1

    # if iter == 119:
    #     filename = f"diffusion_iter{iter:03d}.npz"
    #     load_diffusion = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename, allow_pickle=True)

    #     delta_mom = dif_mom - load_diffusion["mom"]
    #     delta_trac = dif_trac - load_diffusion["trac"]
    #     delta_ke = dif_ke - load_diffusion["ke"]

    #     filename2 = f"velocity_iter{iter:03d}.npz"
    #     load_velocity = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename2, allow_pickle=True)

    #     delta_u_cur = u_cur - load_velocity["u_cur"]
    #     delta_u_bwd = u_bwd - load_velocity["u_bwd"]
    #     delta_u_fwd = u_fwd - load_velocity["u_fwd"]
    #     delta_v_cur = v_cur - load_velocity["v_cur"]
    #     delta_v_bwd = v_bwd - load_velocity["v_bwd"]
    #     delta_v_fwd = v_fwd - load_velocity["v_fwd"]

    #     x=1

    # if iter == 3719:
    #     filename = f"diffusion_iter{iter:04d}.npz"
    #     load_diffusion = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename, allow_pickle=True)

    #     delta_mom = dif_mom - load_diffusion["mom"]
    #     delta_trac = dif_trac - load_diffusion["trac"]
    #     delta_ke = dif_ke - load_diffusion["ke"]

    #     filename2 = f"velocity_iter{iter:04d}.npz"
    #     load_velocity = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + filename2, allow_pickle=True)

    #     delta_u_cur = u_cur - load_velocity["u_cur"]
    #     delta_u_bwd = u_bwd - load_velocity["u_bwd"]
    #     delta_u_fwd = u_fwd - load_velocity["u_fwd"]
    #     delta_v_cur = v_cur - load_velocity["v_cur"]
    #     delta_v_bwd = v_bwd - load_velocity["v_bwd"]
    #     delta_v_fwd = v_fwd - load_velocity["v_fwd"]

    #     x=1

    
    # pom_bgc_1d(iter, base_element, concentration, lambda_w, temp_bwd, sal_bwd, density, ism, swrad, weddy, wgen, wsu, wsv, dz, column_depth, rcp, tracer_map, tracers)
    # pom_bgc_1d(iter, base_element, lambda_w, temp_bwd, sal_bwd, density, ism, swrad, weddy, wgen, wsu, wsv, dif_trac,
    #            dt2, num_layers, z, dz, dzz, dzr, column_depth, 
    #            nrt_o2, nrt_po4, nrt_no3, nrt_nh4, o2b, no3b, ponb_grad, po4b,
    #            smoth, umolbgc, nbcbgc, ntp, rcp, 
    #            concentration, sinking, tracer_map, tracer_type, tracers)

    # Extract concentrations at current and backward time step for leapfrog integration
    conc_cur = concentration[...,iter].copy()
    if iter == 0:   conc_bwd = concentration[...,iter].copy()   # Initialize backward time step at first iteration
    else:           conc_bwd = concentration[...,iter-1].copy()
        
    concentration[...,iter], concentration[...,iter+1] = \
        pom_bgc_1d(iter, configuration, base_element, lambda_w, temp_bwd, sal_bwd, density, ism, swrad, weddy, wgen, wsu, wsv, dif_trac,
                    dt2, num_layers, z, dz, dzz, dzr, column_depth, 
                    nrt_o2, nrt_po4, nrt_no3, nrt_nh4, o2b, no3b, ponb_grad, po4b,
                    smoth, umolbgc, nbcbgc, ntp, rcp, 
                    conc_bwd, conc_cur, sinking, tracer_map, tracer_type, tracers)
    
# ----------------------------------------------------------------------------------------------------
# Write outputs to npz file
# ----------------------------------------------------------------------------------------------------
npp_exists = False  # initialize writing of npp
for trac in tracers:
    if tracers[trac].type == "phytoplankton":
        if not npp_exists:  # first phytoplankton group
            npp = tracers[trac].npp
            npp_exists = True   # npp now exists, update to True to append with npp from later phytoplankton groups
        else:   # subsequent phytoplankton groups
            npp += tracers[trac].npp

conc_daily, conc_monthly = average(concentration,physical,'concentration')
# np.savez('concentration_bfm17-5yr-0907.npz',daily=conc_daily,monthly=conc_monthly)
np.savez('concentration_bfm56-5yr-0907.npz',daily=conc_daily,monthly=conc_monthly)
if npp_exists:
    npp_daily, npp_monthly = average(npp,physical,'npp')
    # np.savez('npp_bfm17-5yr-0907.npz',daily=npp_daily,monthly=npp_monthly)
    np.savez('npp_bfm56-5yr-0907.npz',daily=npp_daily,monthly=npp_monthly)

# np.savez('tracer_indices_bfm17-5yr-0907.npz',**tracer_map)
np.savez('tracer_indices_bfm56-5yr-0907.npz',**tracer_map)


# ----------------------------------------------------------------------------------------------------
# Simulation complete
# ----------------------------------------------------------------------------------------------------
print('Simulation complete.')
elapsed = time.perf_counter() - start
hours, remainder = divmod(elapsed, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"Walltime: {int(hours):02d}:{int(minutes):02d}:{seconds:09.6f}")