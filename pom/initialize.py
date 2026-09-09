import numpy as np
import os
from numba.types import float64, unicode_type
from numba.typed import Dict, List
np.set_printoptions(precision=20)

def initialize_pom(num_layers, temperature_IC, salinity_IC):
    
    # ----------------------------------------------------------------------------------------------------
    # Initialize dictionaries
    # ----------------------------------------------------------------------------------------------------
    # Diffusion ---------------------------------------------------------------
    dif_trac = np.zeros(num_layers,dtype=np.float64)    # Tracers
    dif_mom = np.zeros(num_layers,dtype=np.float64)     # Momentum
    dif_ke = np.zeros(num_layers,dtype=np.float64)      # Kinetic energy

    # Kinetic energy ---------------------------------------------------------------
    ke_cur = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Current
    ke_bwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Backward
    ke_fwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)      # Forward
    # Kinetic energy * length scale
    kel_cur = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Current
    kel_bwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Backward
    kel_fwd = 1.E-07 * np.ones(num_layers,dtype=np.float64)     # Forward

    # Velocity ---------------------------------------------------------------
    # Zonal (u)
    u_cur = np.zeros(num_layers,dtype=np.float64)   # Current
    u_bwd = np.zeros(num_layers,dtype=np.float64)   # Backward
    u_fwd = np.zeros(num_layers,dtype=np.float64)   # Forward
    # Meridional (v)
    v_cur = np.zeros(num_layers,dtype=np.float64)   # Current
    v_bwd = np.zeros(num_layers,dtype=np.float64)   # Backward
    v_fwd = np.zeros(num_layers,dtype=np.float64)   # Forward

    # Temperature ---------------------------------------------------------------
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

    return dif_trac, dif_mom, dif_ke, ke_cur, ke_bwd, ke_fwd, kel_cur, kel_bwd, kel_fwd, \
           u_cur, u_bwd, u_fwd, v_cur, v_bwd, v_fwd, \
           temp_cur, temp_bwd, temp_fwd, temp_adv, temp_int, temp_surf, temp_sflx, temp_bflx, \
           sal_cur, sal_bwd, sal_fwd, sal_adv, sal_int, sal_surf, sal_sflx, sal_bflx, \
           no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb_grad, \
           wsu, wsv, bsu, bsv, ism, swrad, wgen, weddy, mld, \
           counter_ids, counter_params, forcing_ids, forcing_month1, forcing_month2


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

    # forcing_data = [
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
    # ]

    forcing_data = [sclim[:,month], tclim[:,month], wclim[:,month], weddy1[:,month], weddy2[:,month],
                    ism[:,month], wsu[month], wsv[month], swrad[month], wtsurf[month], qcorr[month],
                    no3s[month], nh4s[month], po4s[month], sio4s[month],
                    o2b[month], no3b[month], po4b[month], ponb[month]]

    return forcing_data
