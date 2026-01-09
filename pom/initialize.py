import numpy as np
import os

def initialize_pom(pom1d, physical):
    
    # Initialize dictionaries
    physical["diffusion"] = {
        "tracers": np.zeros(physical["water_column"]["num_layers"]),
        "momentum": np.zeros(physical["water_column"]["num_layers"]),
        "kinetic_energy": np.zeros(physical["water_column"]["num_layers"])
    }
    physical["kinetic_energy"] = {
        # Kinetic energy
        "ke": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
        "keb": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
        "kef": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
        # Kinetic energy * length scale
        "kel": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
        "kelb": 1.E-07 * np.ones(physical["water_column"]["num_layers"]),
        "kelf": 1.E-07 * np.ones(physical["water_column"]["num_layers"])
    }
    physical["velocity"] = {
        # Zonal (U)
        "u": np.zeros(physical["water_column"]["num_layers"]),
        "ub": np.zeros(physical["water_column"]["num_layers"]),
        "uf": np.zeros(physical["water_column"]["num_layers"]),
        # Meridional (V)
        "v": np.zeros(physical["water_column"]["num_layers"]),
        "vb": np.zeros(physical["water_column"]["num_layers"]),
        "vf": np.zeros(physical["water_column"]["num_layers"])
    }
    physical["temperature"] = {
        "t": np.zeros(physical["water_column"]["num_layers"]),      # Current
        "tb": np.zeros(physical["water_column"]["num_layers"]),     # Backward
        "tf": np.zeros(physical["water_column"]["num_layers"]),     # Forward
        "ti": np.zeros(physical["water_column"]["num_layers"]),     # Interpolated (prognostic mode)
        "adv": np.zeros(physical["water_column"]["num_layers"]),    # Lateral advection
        "surf": 0.,         # Surface value
        "surf_flux": 0.,    # Surface flux
        "bot_flux": 0.,     # Bottom flux
    }
    physical["salinity"] = {
        "s": np.zeros(physical["water_column"]["num_layers"]),      # Current
        "sb": np.zeros(physical["water_column"]["num_layers"]),     # Backward
        "sf": np.zeros(physical["water_column"]["num_layers"]),     # Forward
        "si": np.zeros(physical["water_column"]["num_layers"]),     # Interpolated (prognostic mode)
        "adv": np.zeros(physical["water_column"]["num_layers"]),    # Lateral advection
        "surf": 0.,         # Surface value
        "surf_flux": 0.,    # Surface flux
        "bot_flux": 0.,     # Bottom flux
    }
    physical["nutrients"] = {
        # Surface
        "no3s": 0.,         # Nitrate
        "nh4s": 0.,         # Ammonium
        "po4s": 0.,         # Phosphate
        "sio4s": 0.,        # Silicate
        # Bottom
        "o2b": 0.,          # Oxygen
        "no3b": 0.,         # Nitrate
        "po4b": 0.,         # Phosphate
        "ponb_grad": 0.     # PON gradient
    }
    physical["stresses"] = {
        # Wind
        "wsu": 0.,           # Zonal
        "wsv": 0.,           # Meridional
        # Bottom boundary layer
        "bsu": 0.,           # Zonal
        "bsv": 0.            # Meridional
    }
    physical["ism"] = np.zeros(physical["water_column"]["num_layers"])  # Inorganic suspended matter
    physical["swrad"] = 0.      # Shortwave radiation
    physical["wgen"] = np.zeros(physical["water_column"]["num_layers"])
    physical["weddy"] = np.zeros(physical["water_column"]["num_layers"])
    physical["mld"] = np.zeros(physical["water_column"]["num_layers"])  # Mixed layer depth
    
    # Read initial conditions
    # Temperature
    physical["temperature"]["t"] = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])
    physical["temperature"]["tb"] = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])

    # Salinity
    physical["salinity"]["s"] = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])
    physical["salinity"]["sb"] = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])

    # ----------------------------------------------------------------------------------------------------
    # Initialize forcing data
    # ----------------------------------------------------------------------------------------------------
    forcing = {
        "counters": {
            "day_counter": 0,
            "day_interpolator": 0,
            "day_ratio": 0,
            "month_counter": 0,
            "month_interpolator": 0,
            "month_ratio": 0,
            "timesteps_per_day": 0,
            "timesteps_per_month": 0
        },
        "month1": {
            "sclim": np.zeros(physical["water_column"]["num_layers"]),  # Salinity climatology
            "tclim": np.zeros(physical["water_column"]["num_layers"]),  # Temperature climatology
            "wclim": np.zeros(physical["water_column"]["num_layers"]),  # W velocity climatology
            "weddy1": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 1
            "weddy2": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 2
            "ism": np.zeros(physical["water_column"]["num_layers"]),    # Inorganic suspended matter
            "wsu": 0,       # Zonal (U) velocity
            "wsv": 0,       # Meridional (V) velocity
            "swrad": 0,     # Shortwave radiation
            "wtsurf": 0,    # Surface heat flux
            "qcorr": 0,     # Kinetic energy loss
            "no3s": 0,      # Surface no3
            "nh4s": 0,      # Surface nh4
            "po4s": 0,      # Surface po4
            "sio4s": 0,     # Surface sio4
            "o2b": 0,       # Bottom o2
            "no3b": 0,      # Bottom no3
            "po4b": 0,      # Bottom po4
            "ponb": 0       # Bottom pon
        },
        "month2": {
            "sclim": np.zeros(physical["water_column"]["num_layers"]),  # Salinity climatology
            "tclim": np.zeros(physical["water_column"]["num_layers"]),  # Temperature climatology
            "wclim": np.zeros(physical["water_column"]["num_layers"]),  # W velocity climatology
            "weddy1": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 1
            "weddy2": np.zeros(physical["water_column"]["num_layers"]), # Intermittant eddy w velocity 2
            "ism": np.zeros(physical["water_column"]["num_layers"]),    # Inorganic suspended matter
            "wsu": 0,       # Zonal (U) velocity
            "wsv": 0,       # Meridional (V) velocity
            "swrad": 0,     # Shortwave radiation
            "wtsurf": 0,    # Surface heat flux
            "qcorr": 0,     # Kinetic energy loss
            "no3s": 0,      # Surface no3
            "nh4s": 0,      # Surface nh4
            "po4s": 0,      # Surface po4
            "sio4s": 0,     # Surface sio4
            "o2b": 0,       # Bottom o2
            "no3b": 0,      # Bottom no3
            "po4b": 0,      # Bottom po4
            "ponb": 0       # Bottom pon
        }
    }


    return physical, forcing
    

def read_pom_inputs(physical, pom1d):
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
    ism   = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            ism[x,i] = inorganic_suspended_matter_data[physical["water_column"]["num_layers"] * i + x]

    # Surface salinity
    surface_salinity = np.fromfile(os.getcwd() + pom1d["input_files"]["surface_salinity"])

    # Salinity climatology (diagnostic mode)
    salinity_vertical_profile_data = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity"])
    sclim = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            sclim[x,i] = salinity_vertical_profile_data[physical["water_column"]["num_layers"] * i + x]

    # Salinity IC
    salinity = np.fromfile(os.getcwd() + pom1d["input_files"]["salinity_IC"])

    # Temperature climatology (diagnostic mode)
    temperature_vertical_profile_data = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature"])
    tclim = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            tclim[x,i] = temperature_vertical_profile_data[physical["water_column"]["num_layers"] * i + x]

    # Temperature IC
    temperature = np.fromfile(os.getcwd() + pom1d["input_files"]["temperature_IC"])

    # General circulation w velocity climatology
    general_circulation_w_velocity_data = np.fromfile(os.getcwd() + pom1d["input_files"]["w_velocity"])
    wclim  = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            wclim[x,i] = general_circulation_w_velocity_data[physical["water_column"]["num_layers"] * i + x]

    # Intermittant eddy w velocity 1
    intermediate_eddy_w_velocity_1_data = np.fromfile(os.getcwd() + pom1d["input_files"]["eddy_w_velocity_1"])
    weddy1  = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            weddy1[x,i] = intermediate_eddy_w_velocity_1_data[physical["water_column"]["num_layers"] * i + x]


    # Intermittant eddy w velocity 2
    intermediate_eddy_w_velocity_2_data = np.fromfile(os.getcwd() + pom1d["input_files"]["eddy_w_velocity_2"])
    weddy2  = np.zeros((physical["water_column"]["num_layers"],array_length))
    for i in range(0,array_length):
        for x in range(0, physical["water_column"]["num_layers"]):
            weddy2[x,i] = intermediate_eddy_w_velocity_2_data[physical["water_column"]["num_layers"] * i + x]

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

    forcing_data = {
        "sclim": sclim,         # Salinity climatology
        "tclim": tclim,         # Temperature climatology
        "wclim": wclim,         # W velocity climatology
        "weddy1": weddy1,       # Intermittant eddy w velocity 1
        "weddy2": weddy2,       # Intermittant eddy w velocity 1
        "ism": ism,             # Inorganic suspended matter
        "wsu": wsu,             # Zonal (U) velocity
        "wsv": wsv,             # Meridional (V) velocity
        "swrad": swrad,         # Shortwave radiation
        "wtsurf": wtsurf,       # Surface heat flux
        "qcorr": qcorr,         # Kinetic energy loss
        "no3s": no3s,           # Surface no3
        "nh4s": nh4s,           # Surface nh4
        "po4s": po4s,           # Surface po4
        "sio4s": sio4s,         # Surface sio4
        "o2b": o2b,             # Bottom o2
        "no3b": no3b,           # Bottom no3
        "po4b": po4b,           # Bottom po4
        "ponb": ponb            # Bottom pon
    }

    return forcing_data

    
