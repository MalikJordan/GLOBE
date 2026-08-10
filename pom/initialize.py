import numpy as np
import os
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
np.set_printoptions(precision=20)

# @njit
def initialize_pom(num_layers, temperature_IC, salinity_IC):
    
    # ----------------------------------------------------------------------------------------------------
    # Initialize dictionaries
    # ----------------------------------------------------------------------------------------------------
    # Diffusion ---------------------------------------------------------------
    # diffusion = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # diffusion["tracers"] = np.zeros(num_layers,dtype=np.float64)
    # diffusion["momentum"] = np.zeros(num_layers,dtype=np.float64)
    # diffusion["kinetic_energy"] = np.zeros(num_layers,dtype=np.float64)

    dif_trac = np.zeros(num_layers,dtype=np.float64)    # Tracers
    dif_mom = np.zeros(num_layers,dtype=np.float64)     # Momentum
    dif_ke = np.zeros(num_layers,dtype=np.float64)      # Kinetic energy

    # Kinetic energy ---------------------------------------------------------------
    # kinetic_energy = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # # Kinetic energy
    # kinetic_energy["ke"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)
    # kinetic_energy["keb"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)
    # kinetic_energy["kef"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)
    # # Kinetic energy * length scale
    # kinetic_energy["kel"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)
    # kinetic_energy["kelb"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)
    # kinetic_energy["kelf"] = 1.E-07 * np.ones(num_layers,dtype=np.float64)

    # Kinetic energy
    ke_cur = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Current
    ke_bwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Backward
    ke_fwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Forward
    # Kinetic energy * length scale
    kel_cur = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Current
    kel_bwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Backward
    kel_fwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Forward

    # Velocity ---------------------------------------------------------------
    # velocity = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # # Zonal (u)
    # velocity["u"] = np.zeros(num_layers,dtype=np.float64)
    # velocity["ub"] = np.zeros(num_layers,dtype=np.float64)
    # velocity["uf"] = np.zeros(num_layers,dtype=np.float64)
    # # Meridional (v)
    # velocity["v"] = np.zeros(num_layers,dtype=np.float64)
    # velocity["vb"] = np.zeros(num_layers,dtype=np.float64)
    # velocity["vf"] = np.zeros(num_layers,dtype=np.float64)

    # Zonal (u)
    u_cur = np.zeros(num_layers,dtype=np.float64)   # Current
    u_bwd = np.zeros(num_layers,dtype=np.float64)   # Backward
    u_fwd = np.zeros(num_layers,dtype=np.float64)   # Forward
    # Meridional (v)
    v_cur = np.zeros(num_layers,dtype=np.float64)   # Current
    v_bwd = np.zeros(num_layers,dtype=np.float64)   # Backward
    v_fwd = np.zeros(num_layers,dtype=np.float64)   # Forward

    # Temperature ---------------------------------------------------------------
    # temperature = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # temperature["t"] = np.zeros(num_layers,dtype=np.float64)        # Current
    # temperature["tb"] = np.zeros(num_layers,dtype=np.float64)       # Backward
    # temperature["tf"] = np.zeros(num_layers,dtype=np.float64)       # Forward
    # temperature["ti"] = np.zeros(num_layers,dtype=np.float64)       # Interpolated (prognostic mode)
    # temperature["adv"] = np.zeros(num_layers,dtype=np.float64)      # Lateral Advection
    # temperature["surf"] = np.zeros(1,dtype=np.float64)              # Surface value
    # temperature["surf_flux"] = np.zeros(1,dtype=np.float64)         # Surface flux
    # temperature["bot_flux"] = np.zeros(1,dtype=np.float64)          # Bottom flux

    # # Read initial conditions
    # temperature["t"] = np.fromfile(os.getcwd() + temperature_IC)
    # temperature["tb"] = np.fromfile(os.getcwd() + temperature_IC)

    temp_cur = np.zeros(num_layers,dtype=np.float64)                # Current
    temp_bwd = np.zeros(num_layers,dtype=np.float64)                # Backward
    temp_fwd = np.zeros(num_layers,dtype=np.float64)                # Forward
    temp_int = np.zeros(num_layers,dtype=np.float64)                # Interpolated (prognostic mode)
    temp_adv = np.zeros(num_layers,dtype=np.float64)                # Lateral Advection
    temp_surf = np.float64(0.)                                      # Surface value
    temp_sflx = np.float64(0.)                                      # Surface flux
    temp_bflx = np.float64(0.)                                      # Bottom flux

    # Read initial conditions
    temp_cur = np.fromfile(os.getcwd() + temperature_IC)
    temp_bwd = np.fromfile(os.getcwd() + temperature_IC)

    # Salinity ---------------------------------------------------------------
    # salinity = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # salinity["s"] = np.zeros(num_layers,dtype=np.float64)           # Current
    # salinity["sb"] = np.zeros(num_layers,dtype=np.float64)          # Backward
    # salinity["sf"] = np.zeros(num_layers,dtype=np.float64)          # Forward
    # salinity["si"] = np.zeros(num_layers,dtype=np.float64)          # Interpolated (prognostic mode)
    # salinity["adv"] = np.zeros(num_layers,dtype=np.float64)         # Lateral Advection
    # salinity["surf"] = np.zeros(1,dtype=np.float64)                 # Surface value
    # salinity["surf_flux"] = np.zeros(1,dtype=np.float64)            # Surface flux
    # salinity["bot_flux"] = np.zeros(1,dtype=np.float64)             # Bottom flux

    # # Read initial conditions
    # salinity["s"] = np.fromfile(os.getcwd() + salinity_IC)
    # salinity["sb"] = np.fromfile(os.getcwd() + salinity_IC)

    sal_cur = np.zeros(num_layers,dtype=np.float64)                # Current
    sal_bwd = np.zeros(num_layers,dtype=np.float64)                # Backward
    sal_fwd = np.zeros(num_layers,dtype=np.float64)                # Forward
    sal_int = np.zeros(num_layers,dtype=np.float64)                # Interpolated (prognostic mode)
    sal_adv = np.zeros(num_layers,dtype=np.float64)                # Lateral Advection
    sal_surf = np.float64(0.)                                      # Surface value
    sal_sflx = np.float64(0.)                                      # Surface flux
    sal_bflx = np.float64(0.)                                      # Bottom flux

    # Read initial conditions
    sal_cur = np.fromfile(os.getcwd() + salinity_IC)
    sal_bwd = np.fromfile(os.getcwd() + salinity_IC)
    
    # Nutrients ---------------------------------------------------------------
    # nutrients = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    # # Surface
    # nutrients["no3s"] = np.float64(0.)          # Nitrate
    # nutrients["nh4s"] = np.float64(0.)          # Ammonium
    # nutrients["po4s"] = np.float64(0.)          # Phosphate
    # nutrients["sio4s"] = np.float64(0.)         # Silicate
    # # Bottom
    # nutrients["o2b"] = np.float64(0.)           # Oxygen
    # nutrients["no3b"] = np.float64(0.)          # Nitrate
    # nutrients["po4b"] = np.float64(0.)          # Phosphate
    # nutrients["ponb_grad"] = np.float64(0.)     # PON gradient

    # Surface
    no3s = np.float64(0.)                       # Nitrate
    nh4s = np.float64(0.)                       # Ammonium
    po4s = np.float64(0.)                       # Phosphate
    sio4s = np.float64(0.)                      # Silicate

    # Bottom
    o2b = np.float64(0.)                        # Oxygen
    no3b = np.float64(0.)                       # Nitrate
    po4b = np.float64(0.)                       # Phosphate
    ponb_grad = np.float64(0.)                  # PON gradient
    
    # Stresses ---------------------------------------------------------------
    # stresses = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    # # Wind
    # stresses["wsu"] = np.float64(0.)    # Zonal
    # stresses["wsv"] = np.float64(0.)    # Meridional
    # # Bottom boundary layer
    # stresses["bsu"] = np.float64(0.)    # Zonal
    # stresses["bsv"] = np.float64(0.)    # Meridional

    # Wind
    wsu = np.float64(0.)    # Zonal
    wsv = np.float64(0.)    # Meridional
    # Bottom boundary layer
    bsu = np.float64(0.)    # Zonal
    bsv = np.float64(0.)    # Meridional

    # Other ---------------------------------------------------------------
    ism = np.zeros(num_layers, dtype=np.float64)    # Inorganic suspended matter
    swrad = np.float64(0.)                          # Shortwave radiation
    wgen = np.zeros(num_layers, dtype=np.float64)
    weddy = np.zeros(num_layers, dtype=np.float64)
    mld = np.zeros(num_layers, dtype=np.float64)    # Mixed layer depth
    
    # ----------------------------------------------------------------------------------------------------
    # Initialize forcing data
    # ----------------------------------------------------------------------------------------------------
    # counters = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    # counters["day_counter"] = np.float64(0.)
    # counters["day_interpolator"] = np.float64(0.)
    # counters["day_ratio"] = np.float64(0.)
    # counters["month_counter"] = np.float64(0.)
    # counters["month_interpolator"] = np.float64(0.)
    # counters["month_ratio"] = np.float64(0.)
    # counters["timesteps_per_day"] = np.float64(0.)
    # counters["timesteps_per_month"] = np.float64(0.)

    # counters[0] = day_counter
    # counters[1] = day_interpolator
    # counters[2] = day_ratio, 
    # counters[3] = month_counter
    # counters[4] = month_interpolator
    # counters[5] = month_ratio
    # counters[6] = timesteps_per_day
    # counters[7] = timesteps_per_month

    counter_ids = ["day_counter", "day_interpolator", "day_ratio", "month_counter", "month_interpolator", "month_ratio", "timesteps_per_day", "timesteps_per_month"]
    counter_params = [np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.)]

    forcing_ids = ["sclim", "tclim", "wclim", "weddy1", "weddy2", "ism", "wsu", "wsv", "swrad", "wtsurf", "qcorr", "no3s", "nh4s", "po4s", "sio4s", "o2b", "no3b", "po4b", "ponb"]
    forcing_month1 = [np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64),
                      np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.)]
    forcing_month2 = [np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64), np.zeros(num_layers, dtype=np.float64),
                      np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.), np.float64(0.)]

    # forcing[0] = sclim    (Salinity climatology)
    # forcing[1] = tclim    (Temperature climatology)
    # forcing[2] = wclim    (W velocity climatology)
    # forcing[3] = weddy1   (Intermittant eddy w velocity 1)
    # forcing[4] = weddy2   (Intermittant eddy w velocity 2)
    # forcing[5] = ism      (Inorganic suspended matter)
    # forcing[6] = wsu      (Zonal (U) velocity)
    # forcing[7] = wsv      (Meridional (V) velocity)
    # forcing[8] = swrad    (Shortwave radiation)
    # forcing[9] = wtsurf   (Surface heat flux)
    # forcing[10] = qcorr   (Kinetic energy loss)
    # forcing[11] = no3s    (Surface no3)
    # forcing[12] = nh4s    (Surface nh4)
    # forcing[13] = po4s    (Surface po4)
    # forcing[14] = sio4s   (Surface sio4)
    # forcing[15] = o2b     (Bottom o2)
    # forcing[16] = no3b    (Bottom no3)
    # forcing[17] = po4b    (Bottom po4)
    # forcing[18] = ponb    (Bottom pon)

    # month1_arrays = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # month1_arrays["sclim"] = np.zeros(num_layers, dtype=np.float64)     # Salinity climatology
    # month1_arrays["tclim"] = np.zeros(num_layers, dtype=np.float64)     # Temperature climatology
    # month1_arrays["wclim"] = np.zeros(num_layers, dtype=np.float64)     # W velocity climatology
    # month1_arrays["weddy1"] = np.zeros(num_layers, dtype=np.float64)    # Intermittant eddy w velocity 1
    # month1_arrays["weddy2"] = np.zeros(num_layers, dtype=np.float64)    # Intermittant eddy w velocity 2
    # month1_arrays["ism"] = np.zeros(num_layers, dtype=np.float64)       # Inorganic suspended matter

    # month1_floats = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    # month1_floats["wsu"] = np.float64(0.)       # Zonal (U) velocity
    # month1_floats["wsv"] = np.float64(0.)       # Meridional (V) velocity
    # month1_floats["swrad"] = np.float64(0.)     # Shortwave radiation
    # month1_floats["wtsurf"] = np.float64(0.)    # Surface heat flux
    # month1_floats["qcorr"] = np.float64(0.)     # Kinetic energy loss
    # month1_floats["no3s"] = np.float64(0.)      # Surface no3
    # month1_floats["nh4s"] = np.float64(0.)      # Surface nh4
    # month1_floats["po4s"] = np.float64(0.)      # Surface po4
    # month1_floats["sio4s"] = np.float64(0.)     # Surface sio4
    # month1_floats["o2b"] = np.float64(0.)       # Bottom o2
    # month1_floats["no3b"] = np.float64(0.)      # Bottom no3
    # month1_floats["po4b"] = np.float64(0.)      # Bottom po4
    # month1_floats["ponb"] = np.float64(0.)      # Bottom pon

    # month2_arrays = Dict.empty(key_type=types.unicode_type, value_type=types.float64[:])
    # month2_arrays["sclim"] = np.zeros(num_layers, dtype=np.float64)     # Salinity climatology
    # month2_arrays["tclim"] = np.zeros(num_layers, dtype=np.float64)     # Temperature climatology
    # month2_arrays["wclim"] = np.zeros(num_layers, dtype=np.float64)     # W velocity climatology
    # month2_arrays["weddy1"] = np.zeros(num_layers, dtype=np.float64)    # Intermittant eddy w velocity 1
    # month2_arrays["weddy2"] = np.zeros(num_layers, dtype=np.float64)    # Intermittant eddy w velocity 2
    # month2_arrays["ism"] = np.zeros(num_layers, dtype=np.float64)       # Inorganic suspended matter

    # month2_floats = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    # month2_floats["wsu"] = np.float64(0.)       # Zonal (U) velocity
    # month2_floats["wsv"] = np.float64(0.)       # Meridional (V) velocity
    # month2_floats["swrad"] = np.float64(0.)     # Shortwave radiation
    # month2_floats["wtsurf"] = np.float64(0.)    # Surface heat flux
    # month2_floats["qcorr"] = np.float64(0.)     # Kinetic energy loss
    # month2_floats["no3s"] = np.float64(0.)      # Surface no3
    # month2_floats["nh4s"] = np.float64(0.)      # Surface nh4
    # month2_floats["po4s"] = np.float64(0.)      # Surface po4
    # month2_floats["sio4s"] = np.float64(0.)     # Surface sio4
    # month2_floats["o2b"] = np.float64(0.)       # Bottom o2
    # month2_floats["no3b"] = np.float64(0.)      # Bottom no3
    # month2_floats["po4b"] = np.float64(0.)      # Bottom po4
    # month2_floats["ponb"] = np.float64(0.)      # Bottom pon

    # return diffusion, kinetic_energy, velocity, temperature, salinity, nutrients, stresses, ism, swrad, wgen, weddy, mld, \
    #        counters, month1_arrays, month1_floats, month2_arrays, month2_floats
    return dif_trac, dif_mom, dif_ke, ke_cur, ke_bwd, ke_fwd, kel_cur, kel_bwd, kel_fwd, \
           u_cur, u_bwd, u_fwd, v_cur, v_bwd, v_fwd, \
           temp_cur, temp_bwd, temp_fwd, temp_adv, temp_int, temp_surf, temp_sflx, temp_bflx, \
           sal_cur, sal_bwd, sal_fwd, sal_adv, sal_int, sal_surf, sal_sflx, sal_bflx, \
           no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb_grad, \
           wsu, wsv, bsu, bsv, ism, swrad, wgen, weddy, mld, \
           counter_ids, counter_params, forcing_ids, forcing_month1, forcing_month2


# def initialize_pom(pom1d, physical):
    
#     # Initialize dictionaries
#     physical["diffusion"] = {
#         "tracers": np.zeros(physical["water_column"]["num_layers"]),
#         "momentum": np.zeros(physical["water_column"]["num_layers"]),
#         "kinetic_energy": np.zeros(physical["water_column"]["num_layers"])
#     }
#     physical["kinetic_energy"] = {
#         # Kinetic energy
#         "ke": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
#         "keb": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
#         "kef": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
#         # Kinetic energy * length scale
#         "kel": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
#         "kelb": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
#         "kelf": 1.E-07 * np.ones(physical["water_column"]["num_layers"])
#     }
#     physical["velocity"] = {
#         # Zonal (U)
#         "u": np.zeros(physical["water_column"]["num_layers"]),
#         "ub": np.zeros(physical["water_column"]["num_layers"]),
#         "uf": np.zeros(physical["water_column"]["num_layers"]),
#         # Meridional (V)
#         "v": np.zeros(physical["water_column"]["num_layers"]),
#         "vb": np.zeros(physical["water_column"]["num_layers"]),
#         "vf": np.zeros(physical["water_column"]["num_layers"])
#     }
#     physical["temperature"] = {
#         "t": np.zeros(physical["water_column"]["num_layers"]),      # Current
#         "tb": np.zeros(physical["water_column"]["num_layers"]),     # Backward
#         "tf": np.zeros(physical["water_column"]["num_layers"]),     # Forward
#         "ti": np.zeros(physical["water_column"]["num_layers"]),     # Interpolated (prognostic mode)
#         "adv": np.zeros(physical["water_column"]["num_layers"]),    # Lateral advection
#         "surf": 0.,         # Surface value
#         "surf_flux": 0.,    # Surface flux
#         "bot_flux": 0.,     # Bottom flux
#     }
#     physical["salinity"] = {
#         "s": np.zeros(physical["water_column"]["num_layers"]),      # Current
#         "sb": np.zeros(physical["water_column"]["num_layers"]),     # Backward
#         "sf": np.zeros(physical["water_column"]["num_layers"]),     # Forward
#         "si": np.zeros(physical["water_column"]["num_layers"]),     # Interpolated (prognostic mode)
#         "adv": np.zeros(physical["water_column"]["num_layers"]),    # Lateral advection
#         "surf": 0.,         # Surface value
#         "surf_flux": 0.,    # Surface flux
#         "bot_flux": 0.,     # Bottom flux
#     }
#     physical["nutrients"] = {
#         # Surface
#         "no3s": 0.,         # Nitrate
#         "nh4s": 0.,         # Ammonium
#         "po4s": 0.,         # Phosphate
#         "sio4s": 0.,        # Silicate
#         # Bottom
#         "o2b": 0.,          # Oxygen
#         "no3b": 0.,         # Nitrate
#         "po4b": 0.,         # Phosphate
#         "ponb_grad": 0.     # PON gradient
#     }
#     physical["stresses"] = {
#         # Wind
#         "wsu": 0.,           # Zonal
#         "wsv": 0.,           # Meridional
#         # Bottom boundary layer
#         "bsu": 0.,           # Zonal
#         "bsv": 0.            # Meridional
#     }
#     physical["ism"] = np.zeros(physical["water_column"]["num_layers"])  # Inorganic suspended matter
#     physical["swrad"] = 0.      # Shortwave radiation
#     physical["wgen"] = np.zeros(physical["water_column"]["num_layers"])
#     physical["weddy"] = np.zeros(physical["water_column"]["num_layers"])
#     physical["mld"] = np.zeros(physical["water_column"]["num_layers"])  # Mixed layer depth
    
#     # Read initial conditions
#     # Temperature
#     physical["temperature"]["t"] = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])
#     physical["temperature"]["tb"] = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])

#     # Salinity
#     physical["salinity"]["s"] = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])
#     physical["salinity"]["sb"] = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])

#     # ----------------------------------------------------------------------------------------------------
#     # Initialize forcing data
#     # ----------------------------------------------------------------------------------------------------
#     forcing = {
#         "counters": {
#             "day_counter": 0,
#             "day_interpolator": 0,
#             "day_ratio": 0,
#             "month_counter": 0,
#             "month_interpolator": 0,
#             "month_ratio": 0,
#             "timesteps_per_day": 0,
#             "timesteps_per_month": 0
#         },
#         "month1": {
#             "sclim": np.zeros(physical["water_column"]["num_layers"]),  # Salinity climatology
#             "tclim": np.zeros(physical["water_column"]["num_layers"]),  # Temperature climatology
#             "wclim": np.zeros(physical["water_column"]["num_layers"]),  # W velocity climatology
#             "weddy1": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 1
#             "weddy2": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 2
#             "ism": np.zeros(physical["water_column"]["num_layers"]),    # Inorganic suspended matter
#             "wsu": 0,       # Zonal (U) velocity
#             "wsv": 0,       # Meridional (V) velocity
#             "swrad": 0,     # Shortwave radiation
#             "wtsurf": 0,    # Surface heat flux
#             "qcorr": 0,     # Kinetic energy loss
#             "no3s": 0,      # Surface no3
#             "nh4s": 0,      # Surface nh4
#             "po4s": 0,      # Surface po4
#             "sio4s": 0,     # Surface sio4
#             "o2b": 0,       # Bottom o2
#             "no3b": 0,      # Bottom no3
#             "po4b": 0,      # Bottom po4
#             "ponb": 0       # Bottom pon
#         },
#         "month2": {
#             "sclim": np.zeros(physical["water_column"]["num_layers"]),  # Salinity climatology
#             "tclim": np.zeros(physical["water_column"]["num_layers"]),  # Temperature climatology
#             "wclim": np.zeros(physical["water_column"]["num_layers"]),  # W velocity climatology
#             "weddy1": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 1
#             "weddy2": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 2
#             "ism": np.zeros(physical["water_column"]["num_layers"]),    # Inorganic suspended matter
#             "wsu": 0,       # Zonal (U) velocity
#             "wsv": 0,       # Meridional (V) velocity
#             "swrad": 0,     # Shortwave radiation
#             "wtsurf": 0,    # Surface heat flux
#             "qcorr": 0,     # Kinetic energy loss
#             "no3s": 0,      # Surface no3
#             "nh4s": 0,      # Surface nh4
#             "po4s": 0,      # Surface po4
#             "sio4s": 0,     # Surface sio4
#             "o2b": 0,       # Bottom o2
#             "no3b": 0,      # Bottom no3
#             "po4b": 0,      # Bottom po4
#             "ponb": 0       # Bottom pon
#         }
#     }


#     return physical, forcing

    

def read_pom_inputs(num_layers, pom1d, month):
    """
    Description: Opens forcing files reading the paths specified in the pom_input namelist.

    :return: data arrays for wind stress, surface salinity, solar radiation, inorganic
             suspended matter, salinity and temperature vertical profiles, general circulation
             for w velocity, intermediate eddy velocities, salinity and temperature initial
             conditions, heat flux loss, and surface and bottom nutrients
    """

    # Length of input arrays
    array_length = 13   # months (D-J-F-M-A-M-J-J-A-S-O-N-D)
    
    # Wind speed
    wind_speed_data = np.fromfile(os.getcwd() + pom1d["input_files"]["wind_stress"])
    wsu   = np.zeros(array_length)
    wsv   = np.zeros(array_length)
    for i in range(0,array_length):
        wsu[i] = wind_speed_data[2*i + 0]
        wsv[i] = wind_speed_data[2*i + 1]

    # SOlar radiation, shortwave radiation, and heat flux
    solar_radiation = np.fromfile(os.getcwd() + pom1d["input_files"]["shortwave_solar_radiation"])
    heat_flux_loss_data = np.fromfile(os.getcwd() + pom1d["input_files"]["heat_flux_loss"])
    swrad = np.zeros(array_length)
    wtsurf = np.zeros(array_length)
    qcorr = np.zeros(array_length)
    for i in range(0,array_length):
        swrad[i]  = heat_flux_loss_data[3*i + 0]
        wtsurf[i] = heat_flux_loss_data[3*i + 1]
        qcorr[i]  = heat_flux_loss_data[3*i + 2]

    # Inorganic suspended matter
    inorganic_suspended_matter_data = np.fromfile(os.getcwd() + pom1d["input_files"]["inorganic_suspended_matter"])
    ism   = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            ism[x,i] = inorganic_suspended_matter_data[num_layers * i + x]

    # Surface salinity
    surface_salinity = np.fromfile(os.getcwd() + pom1d["input_files"]["surface_salinity"])

    # Salinity climatology (diagnostic mode)
    salinity_vertical_profile_data = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity"])
    sclim = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            sclim[x,i] = salinity_vertical_profile_data[num_layers * i + x]

    # Salinity IC
    salinity = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])

    # Temperature climatology (diagnostic mode)
    temperature_vertical_profile_data = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature"])
    tclim = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            tclim[x,i] = temperature_vertical_profile_data[num_layers * i + x]

    # Temperature IC
    temperature = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])

    # General circulation w velocity climatology
    general_circulation_w_velocity_data = np.fromfile(os.getcwd() + pom1d["input_files"]["w_velocity"])
    wclim  = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            wclim[x,i] = general_circulation_w_velocity_data[num_layers * i + x]

    # Intermittant eddy w velocity 1
    intermediate_eddy_w_velocity_1_data = np.fromfile(os.getcwd() + pom1d["input_files"]["eddy_w_velocity_1"])
    weddy1  = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            weddy1[x,i] = intermediate_eddy_w_velocity_1_data[num_layers * i + x]


    # Intermittant eddy w velocity 2
    intermediate_eddy_w_velocity_2_data = np.fromfile(os.getcwd() + pom1d["input_files"]["eddy_w_velocity_2"])
    weddy2  = np.zeros((num_layers,array_length))
    for i in range(0,array_length):
        for x in range(0, num_layers):
            weddy2[x,i] = intermediate_eddy_w_velocity_2_data[num_layers * i + x]

    # Surface nutrients
    surface_nutrients_data  = np.fromfile(os.getcwd() + pom1d["input_files"]["surface_nutrients"])
    no3s  = np.zeros(array_length)
    nh4s  = np.zeros(array_length)
    po4s  = np.zeros(array_length)
    sio4s = np.zeros(array_length)
    for i in range(0,array_length):
        no3s[i]  = surface_nutrients_data[4*i + 0]
        nh4s[i]  = surface_nutrients_data[4*i + 1]
        po4s[i]  = surface_nutrients_data[4*i + 2]
        sio4s[i] = surface_nutrients_data[4*i + 3]

    # Bottom nutrients
    bottom_nutrients_data = np.fromfile(os.getcwd() + pom1d["input_files"]["bottom_nutrients"])
    o2b   = np.zeros(array_length)
    no3b  = np.zeros(array_length)
    po4b  = np.zeros(array_length)
    ponb  = np.zeros(array_length)
    for i in range(0,array_length):
        o2b[i]  = bottom_nutrients_data[4*i + 0]
        no3b[i] = bottom_nutrients_data[4*i + 1]
        po4b[i] = bottom_nutrients_data[4*i + 2]
        ponb[i] = bottom_nutrients_data[4*i + 3]

    # forcing_data = {
    #     "sclim": sclim,         # Salinity climatology
    #     "tclim": tclim,         # Temperature climatology
    #     "wclim": wclim,         # W velocity climatology
    #     "weddy1": weddy1,       # Intermittant eddy w velocity 1
    #     "weddy2": weddy2,       # Intermittant eddy w velocity 1
    #     "ism": ism,             # Inorganic suspended matter
    #     "wsu": wsu,             # Zonal (U) velocity
    #     "wsv": wsv,             # Meridional (V) velocity
    #     "swrad": swrad,         # Shortwave radiation
    #     "wtsurf": wtsurf,       # Surface heat flux
    #     "qcorr": qcorr,         # Kinetic energy loss
    #     "no3s": no3s,           # Surface no3
    #     "nh4s": nh4s,           # Surface nh4
    #     "po4s": po4s,           # Surface po4
    #     "sio4s": sio4s,         # Surface sio4
    #     "o2b": o2b,             # Bottom o2
    #     "no3b": no3b,           # Bottom no3
    #     "po4b": po4b,           # Bottom po4
    #     "ponb": ponb            # Bottom pon
    # }

    forcing_data = [sclim[:,month], tclim[:,month], wclim[:,month], weddy1[:,month], weddy2[:,month],
                    ism[:,month], wsu[month], wsv[month], swrad[month], wtsurf[month], qcorr[month],
                    no3s[month], nh4s[month], po4s[month], sio4s[month],
                    o2b[month], no3b[month], po4b[month], ponb[month]]

    return forcing_data

    
