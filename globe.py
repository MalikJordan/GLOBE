import os
# import sys
import numpy as np
import yaml
from setup.initialize import import_bgc_model, import_physical_model
from functions.bgc_rate_eqns import bgc_rate_eqns
from pom.calculations import density_profile, kinetic_energy_profile, temperature_and_salinity_profiles, zonal_velocity_profile, meridional_velocity_profile
from pom.forcing import forcing_manager
from pom.initialize import initialize_pom
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

file = 'tests/bfm17/bfm17-1d.yaml'
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
    pom1d["general"]["coriolis"] = 2. * pom1d["general"]["earth_angular_velocity"] * np.sin(physical["environment"]["latitude"] * 2. * np.pi / 360.)

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
        physical["water_column"]["upper_depth"]
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

    # if not pom_bfm_parameters.pom_only:
    #     bfm_phys_vars = pom_to_bfm(bfm_phys_vars, vertical_grid, temperature, salinity, inorganic_suspended_matter, shortwave_radiation, density, wind_stress)
    #     bfm_phys_vars.vertical_extinction = vertical_extinction(bfm_phys_vars, d3state, species)
    #     bfm_phys_vars.irradiance = light_distribution(bfm_phys_vars)

    #     d3state, d3stateb, d3ave = pom_bfm_1d(i, vertical_grid, t, diffusion, nutrients, bfm_phys_vars, d3state, d3stateb, d3ave, include, species)    






concentration = []
tracer_indices = {}
index = 0   # counting number to keep track of tracer index in concentration matrix
for t in tracers:
    num_constituents = len(tracers[t].composition)
    tracer_indices[t] = list(np.arange(index, index+num_constituents, 1))
    for i in range(num_constituents):
        concentration.append(tracers[t].conc[i,...])

        index += 1

concentration = np.array(concentration,dtype=float)

np.savez('npzd.npz',concentration=concentration,time=physical["simulation"]["time"])
np.savez('tracer_indices_npzd.npz',**tracer_indices)

print('Simulation complete.')




month1 = {
    "sclim": np.zeros(physical["water_column"]["num_layers"]),
    "tclim": np.zeros(physical["water_column"]["num_layers"]),
    "wclim": np.zeros(physical["water_column"]["num_layers"]),
    "weddy1": np.zeros(physical["water_column"]["num_layers"]),
    "weddy2": np.zeros(physical["water_column"]["num_layers"]),
    "ism": np.zeros(physical["water_column"]["num_layers"]),
    "wsu": 0,
    "wsv": 0,
    "swrad": 0,
    "wtsurf": 0,
    "qcorr": 0,
    "no3s": 0,
    "nh4s": 0,
    "po4s": 0,
    "sio4s": 0,
    "o2b": 0,
    "no3b": 0,
    "po4b": 0,
    "ponb": 0
}