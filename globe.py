import os
# import sys
import numpy as np
import yaml
from setup.initialize import import_bgc_model, import_physical_model
from functions.bgc_rate_eqns import bgc_rate_eqns
from functions.calculate_averages import average
from pom.calculations import density_profile, kinetic_energy_profile, temperature_and_salinity_profiles, zonal_velocity_profile, meridional_velocity_profile
from pom.forcing import forcing_manager
from pom.initialize import initialize_pom
from pom.coupling import pom_bgc_1d
from pom.check_phys import dens, u, ub, v, vb, t, tb, s, sb, q2, q2b, q2l, q2lb, km, kh, kq
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
file = 'physical.yaml'
file_path = os.getcwd() + '/' + file
physical = import_physical_model(file_path)

file = 'tests/bfm17/bfm17-1d-2.yaml'
file_path = os.getcwd() + '/' + file
base_element, reactions, tracers = import_bgc_model(file_path, physical)

# ----------------------------------------------------------------------------------------------------
# Initialize POM1D (if necessary)
# ----------------------------------------------------------------------------------------------------
if physical["environment"]["forcing"] == "pom1d":
    with open(os.getcwd() + '/pom1d.yaml', 'r') as f:
        pom1d = yaml.full_load(f)
    physical, forcing = initialize_pom(pom1d, physical)
    physical = density_profile(physical)
    # pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * np.pi / 360.)
    pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * (3.14159265359) / 360.)


# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
# for iter in range(0,physical["simulation"]["iters"]-1):
#     bgc_rate_eqns(iter, base_element, parameters, tracers)

for iter in range(0,physical["simulation"]["iters"]-1):
    
    # Turbulence closure
    physical["kinetic_energy"]["kef"][:] = physical["kinetic_energy"]["keb"][:]
    physical["kinetic_energy"]["kelf"][:] = physical["kinetic_energy"]["kelb"][:]
    physical = kinetic_energy_profile(physical, pom1d)

    # Define forcings
    physical, forcing = forcing_manager(iter, physical, forcing, pom1d)

    # Temperature and salinity computation
    if pom1d["general"]["idiagn"] == 0:
        # Prognostic mode
        # Temperature and salinity fully computed by model
        physical["temperature"]["surf"] = physical["temperature"]["tf"][0]
        physical["salinity"]["surf"] = physical["salinity"]["sf"][0]
        if pom1d["relaxation_times"]["trt"] != 0:
            for j in range(0, physical["water_column"]["num_layers"]):
                if (-physical["vertical_grid"]["dzz"][j] * physical["water_column"]["column_depth"]) >= physical["water_column"]["upper_depth"]:
                    physical["temperature"]["adv"][j] = (physical["temperature"]["ti"][j] - physical["temperature"]["t"][j]) / (pom1d["relaxation_times"]["trt"] * physical["simulation"]["sec_per_day"])

        if pom1d["relaxation_times"]["srt"] != 0:
            for j in range(0, physical["water_column"]["num_layers"]):
                if (-physical["vertical_grid"]["dzz"][j] * physical["water_column"]["column_depth"]) >= physical["water_column"]["upper_depth"]:
                    physical["salinity"]["adv"][j] = (physical["salinity"]["si"][j] - physical["salinity"]["s"][j]) / (pom1d["relaxation_times"]["srt"] * physical["simulation"]["sec_per_day"])
        
        # Calculate surface salinity flux
        physical["salinity"]["surf_flux"] = -(physical["salinity"]["surf"] - physical["salinity"]["s"][0]) * pom1d["relaxation_times"]["srt"] / physical["simulation"]["sec_per_day"]

        # Calculate temperature
        physical["temperature"]["tf"][:] = physical["temperature"]["tb"][:] + (physical["temperature"]["adv"][:] * physical["simulation"]["dt2"])
        physical["temperature"] = temperature_and_salinity_profiles(physical, pom1d, physical["temperature"], 'Temperature')

        # Calculate salinity
        physical["salinity"]["sf"][:] = physical["salinity"]["sb"][:] + (physical["salinity"]["adv"][:] * physical["simulation"]["dt2"])
        physical["salinity"] = temperature_and_salinity_profiles(physical, pom1d, physical["salinity"], 'Salinity')

        # Mix the timestep (Asselin filter)
        physical["temperature"]["t"][:] = physical["temperature"]["t"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["temperature"]["tf"] + physical["temperature"]["tb"][:] - 2. * physical["temperature"]["t"])
        physical["salinity"]["s"][:] = physical["salinity"]["s"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["salinity"]["sf"][:] + physical["salinity"]["sb"][:] - 2. * physical["salinity"]["s"][:])

    # Velocity computation
    physical["velocity"]["uf"][:] = physical["velocity"]["ub"][:] + physical["simulation"]["dt2"] * pom1d["general"]["coriolis"] * physical["velocity"]["v"][:]
    physical = zonal_velocity_profile(physical, pom1d)

    physical["velocity"]["vf"][:] = physical["velocity"]["vb"][:] - physical["simulation"]["dt2"] * pom1d["general"]["coriolis"] * physical["velocity"]["u"][:]
    physical = meridional_velocity_profile(physical, pom1d)

    # Mix the timestep (Asselin filter)
    physical["kinetic_energy"]["ke"][:] = physical["kinetic_energy"]["ke"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["kinetic_energy"]["kef"][:] + physical["kinetic_energy"]["keb"][:] - 2. * physical["kinetic_energy"]["ke"][:])
    physical["kinetic_energy"]["kel"][:] = physical["kinetic_energy"]["kel"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["kinetic_energy"]["kelf"][:] + physical["kinetic_energy"]["kelb"][:] - 2. * physical["kinetic_energy"]["kel"][:])

    physical["velocity"]["u"][:] = physical["velocity"]["u"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["velocity"]["uf"][:] + physical["velocity"]["ub"][:] - 2. * physical["velocity"]["u"][:])
    physical["velocity"]["v"][:] = physical["velocity"]["v"][:] + 0.5 * pom1d["general"]["smoth"] * (physical["velocity"]["vf"][:] + physical["velocity"]["vb"][:] - 2. * physical["velocity"]["v"][:])

    # Restore the time sequence
    physical["kinetic_energy"]["keb"][:] = physical["kinetic_energy"]["ke"][:]
    physical["kinetic_energy"]["ke"][:] = physical["kinetic_energy"]["kef"][:]
    physical["kinetic_energy"]["kelb"][:] = physical["kinetic_energy"]["kel"][:]
    physical["kinetic_energy"]["kel"][:] = physical["kinetic_energy"]["kelf"][:]

    physical["velocity"]["ub"][:] = physical["velocity"]["u"][:]
    physical["velocity"]["u"][:] = physical["velocity"]["uf"][:]
    physical["velocity"]["vb"][:] = physical["velocity"]["v"][:]
    physical["velocity"]["v"][:] = physical["velocity"]["vf"][:]

    physical["temperature"]["tb"][:] = physical["temperature"]["t"][:]
    physical["temperature"]["t"][:] = physical["temperature"]["tf"][:]
    physical["salinity"]["sb"][:] = physical["salinity"]["s"][:]
    physical["salinity"]["s"][:] = physical["salinity"]["sf"][:]

    # Update density
    physical = density_profile(physical)

    if iter < 10:
        delta_q2 = physical["kinetic_energy"]["ke"][:150] - q2[iter]
        delta_q2b = physical["kinetic_energy"]["keb"][:150] - q2b[iter]
        delta_q2l = physical["kinetic_energy"]["kel"][:150] - q2l[iter]
        delta_q2lb = physical["kinetic_energy"]["kelb"][:150] - q2lb[iter]

        delta_t = physical["temperature"]["t"][:150] - t[iter]
        delta_tb = physical["temperature"]["tb"][:150] - tb[iter]
        delta_s = physical["salinity"]["s"][:150] - s[iter]
        delta_sb = physical["salinity"]["sb"][:150] - sb[iter]

        delta_u = physical["velocity"]["u"][:150] - u[iter]
        delta_ub = physical["velocity"]["ub"][:150] - ub[iter]
        delta_v = physical["velocity"]["v"][:150] - v[iter]
        delta_vb = physical["velocity"]["vb"][:150] - vb[iter]

        delta_rho = physical["density"][:150] - dens[iter]

        delta_km = physical["diffusion"]["momentum"][:] - km[iter]
        delta_kh = physical["diffusion"]["tracers"][:] - kh[iter]
        delta_kq = physical["diffusion"]["kinetic_energy"][:] - kq[iter]

        x = 1

    pom_bgc_1d(iter, base_element, physical, pom1d, tracers)


# ----------------------------------------------------------------------------------------------------
# Write outputs to npz file
# ----------------------------------------------------------------------------------------------------
concentration = []
npp_exists = False  # initialize writing of npp
tracer_indices = {} # used to keep track of tracer/consitient index in concentration matriz
index = 0   # counting number for tracer indices
for trac in tracers:
    num_constituents = len(tracers[trac].composition)
    tracer_indices[trac] = list(np.arange(index, index+num_constituents, 1))   # identify tracer constituents with their own index
    for i in range(num_constituents):
        concentration.append(tracers[trac].conc[i,...])    # add concentration to matrix

        index += 1  # update index

    # Update npp
    if tracers[trac].type == "phytoplankton":
        if not npp_exists:  # first phytoplankton group
            npp = tracers[trac].npp
            npp_exists = True   # npp now exists, update to True to append with npp from later phytoplankton groups
        else:   # subsequent phytoplankton groups
            npp += tracers[trac].npp

concentration = np.array(concentration,dtype=float) # convert concentration from list to array

conc_daily, conc_monthly = average(concentration,physical,'concentration')
np.savez('concentration.npz',daily=conc_daily,monthly=conc_monthly)
if npp_exists:
    npp_daily, npp_monthly = average(npp,physical,'npp')
    np.savez('npp.npz',daily=npp_daily,monthly=npp_monthly)


# np.savez('bfm17-1d-2-6.npz',concentration=concentration,npp=npp,time=physical["simulation"]["time"])
np.savez('tracer_indices_bfm17-1d-2-6.npz',**tracer_indices)


# ----------------------------------------------------------------------------------------------------
# Simulation complete
# ----------------------------------------------------------------------------------------------------
print('Simulation complete.')
