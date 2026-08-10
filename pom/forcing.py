import numpy as np
import os
from pom.initialize import read_pom_inputs
np.set_printoptions(precision=20)


def forcing_manager(iter, dt, num_layers, pom1d, counters, month1, month2):
    """
    Description: Handles the forcing data reading and interpolation.
                 Accomodates perpetual year monthly forcing time series.
                 Data currently handled are time series of monthly averaged data.
                 
    :return: forcing data
    """
    # counters[0] = day_counter
    # counters[1] = day_interpolator
    # counters[2] = day_ratio, 
    # counters[3] = month_counter
    # counters[4] = month_interpolator
    # counters[5] = month_ratio
    # counters[6] = timesteps_per_day
    # counters[7] = timesteps_per_month

    # forcing[0] = sclim    (Salinity climatology)
    # forcing[1] = tclim    (Temperature climatology)
    # forcing[2] = wclim    (W velocity climatology)
    # forcing[3] = weddy1   (Intermittant eddy w velocity 1)
    # forcing[4] = weddy2   (Intermittant eddy w velocity 2)
    # forcing[5] = ism      (Inorganic suspended matter)
    # forcing[6] = wsu      (Zonal (U) velocity)
    # forcing[7] = wsv      (Meridional (V) velocity)
    # forcing[8] = swrad      (IShortwave radiation)
    # forcing[9] = wtsurf      (Surface heat flux)
    # forcing[10] = qcorr      (Kinetic energy loss)
    # forcing[11] = no3s      (Surface no3)
    # forcing[12] = nh4s      (Surface nh4)
    # forcing[13] = po4s      (Surface po4)
    # forcing[14] = sio4s      (Surface sio4)
    # forcing[15] = o2b      (Bottom o2)
    # forcing[16] = no3b      (Bottom no3)
    # forcing[17] = po4b      (Bottom po4)
    # forcing[18] = ponb      (Bottom pon)
    # INITIALISATION AND FIRST FORCING READING
    if iter == 0:

        # Month and day counters
        counters[0] = 1   # initialize day counter
        counters[3] = 0   # initialize month counter

        # Timesteps to cover one day/month
        counters[6] = 86400. / dt # define time steps for day
        counters[7] = 30 * counters[6]  # define time steps per month
        
        # Month and day interpolators
        counters[1] = -1  # initialize day interpolator
        counters[4] = (counters[7] / 2.) - 1    # initialize month interpolator (climatological forcing data centered at day 15 of each month)

        # Reading for the first month
        month1 = write_forcing_data(counters[3], num_layers, pom1d)

        # Update the day counter
        counters[0] += 1

        # Update the month counter
        counters[3] += 1

        # Reading for the second month
        month2 = write_forcing_data(counters[3], num_layers, pom1d)

    # Update interpolation counters
    counters[1] += 1  # update day interpolator
    counters[2] = counters[1] / counters[6]   # update day ratio

    if counters[2] == 1:    # day ratio == 1
        x=1

    counters[4] += 1  # update month interpolator
    counters[5] = counters[4] / counters[7]   # update month ratio

    # Interpolate wind stress
    wsu = month1[6] + counters[5] * (month2[6] - month1[6]) # zonal wind stress
    wsv = month1[7] + counters[5] * (month2[7] - month1[7]) # meridional wind stress

    # Interpolate heat flux
    if pom1d["general"]["idiagn"] == 0: # prognostic mode (idiagn = 0)
        temp_sflx = month1[9] + counters[5] * (month2[9] - month1[9])   # surface temperature flux
        swrad = month1[8] + counters[5] * (month2[8] - month1[8])       # shortwave radiation
    else:   # diagnostic mode (idiagn = 1)
        temp_sflx = 0.  # surface temperature flux not needed for diagnostic mode
        swrad = month1[8] + counters[5] * (month2[8] - month1[8])       # shortwave radiation

    # Interpolate temperature and salinity profiles
    temp_int = month1[1] + counters[5] * (month2[1] - month1[1])    # interpolated temperature
    sal_int = month1[0] + counters[5] * (month2[0] - month1[0])      # interpolated salinity
    wgen = month1[2] + counters[5] * (month2[2] - month1[2])        # 

    if counters[5] <= 0.5:  # during first 15 days
        weddy = month1[3]   # intermittant eddy w velocity 1
    else:
        weddy = month1[4]   # intermittant eddy w velocity 2

    if pom1d["general"]["idiagn"] == 0:
        temp_surf = temp_int[0] # surface temperature
        sal_surf = sal_int[0]   # surface salinity
    elif pom1d["general"]["idiagn"] == 1:
        temp_surf = 0.
        sal_surf = 0.
        temp_fwd = temp_int # forward temperature profile
        sal_fwd = sal_int   # forwards salinity profile

    # Interpolate suspended inorganic matter
    ism = month1[5] + counters[5] * (month2[5] - month1[5])
    
    # Interpolate surface and bottom nutrients
    no3s = month1[11] + counters[5] * (month2[11] - month1[11])     # surface nitrate
    nh4s = month1[12] + counters[5] * (month2[12] - month1[12])     # surface ammonium
    po4s = month1[13] + counters[5] * (month2[13] - month1[13])     # surface phosphate
    sio4s = month1[14] + counters[5] * (month2[14] - month1[14])    # surface silicate

    o2b = month1[15] + counters[5] * (month2[15] - month1[15])          # bottom oxygen
    no3b = month1[16] + counters[5] * (month2[16] - month1[16])         # bottom nitrate
    po4b = month1[17] + counters[5] * (month2[17] - month1[17])         # bottom phosphate
    ponb_grad = month1[18] + counters[5] * (month2[18] - month1[18])    # bottom pon gradient

    if counters[4] == counters[7]:  # if month interpolator == time steps per month

        # Update the month counter
        if os.path.exists('output_file.txt'):
            output_file = open('output_file.txt','a')
            output_file.write('month_counter = ')
            output_file.write(str(counters[3]))
            output_file.write('\n')
            output_file.close()
        else:
            output_file = open('output_file.txt','w')
            output_file.write('month_counter = ')
            output_file.write(str(counters[3]))
            output_file.write('\n')
            output_file.close()

        print('month_counter = ',counters[3])
        counters[3] += 1    # month counter
        
        # Reset the interpolator
        counters[4] = 0

        # Shift the monthly data
        month1 = month2

        # If 12 months have passed, restart the reading sequence
        if counters[3] > 12:
            counters[3] = 0     # month counter = 0
            month1 = write_forcing_data(counters[3], num_layers, pom1d)
            counters[3] += 1    # update month counter

        # Read the following month
        month2 = write_forcing_data(counters[3], num_layers, pom1d)

    return month1, month2, temp_fwd, temp_int, temp_surf, temp_sflx, sal_fwd, sal_int, sal_surf, \
           ism, swrad, wsu, wsv, wgen, weddy, no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb_grad



def write_forcing_data(month_counter, num_layers, pom1d):
    """
    Definition: Writes forcing data for each month.
    """
    forcing_data = read_pom_inputs(num_layers, pom1d, month_counter)
    
    # Climatology
    # month_data = forcing_data
    # month_data["sclim"] = forcing_data["sclim"][:,month_counter]
    # month_data["tclim"] = forcing_data["tclim"][:,month_counter]
    # month_data["wclim"] = forcing_data["wclim"][:,month_counter]

    # # Intermittant eddy w velocities
    # month_data["weddy1"] = forcing_data["weddy1"][:,month_counter]
    # month_data["weddy2"] = forcing_data["weddy2"][:,month_counter]

    # # Inorganic suspended matter
    # month_data["ism"] = forcing_data["ism"][:-1,month_counter]

    # month_data["wsu"] = -0.001 * forcing_data["wsu"][month_counter]     # N/m2-->m2/s2
    # month_data["wsv"] = -0.001 * forcing_data["wsv"][month_counter]     # N/m2-->m2/s2
    
    # # Heat Flux
    # month_data["swrad"]  = -forcing_data["swrad"][month_counter] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
    # month_data["wtsurf"] = -forcing_data["wtsurf"][month_counter] / pom1d["general"]["water_specific_heat_times_density"]   # W/m2-->deg.C*m/s
    # month_data["qcorr"]  = forcing_data["qcorr"][month_counter]

    # # Surface nutrients
    # month_data["no3s"]  = forcing_data["no3s"][month_counter]
    # month_data["nh4s"]  = forcing_data["nh4s"][month_counter]
    # month_data["po4s"]  = forcing_data["po4s"][month_counter]
    # month_data["sio4s"] = forcing_data["sio4s"][month_counter]

    # # Bottom nutrients
    # month_data["o2b"]  = forcing_data["o2b"][month_counter]
    # month_data["no3b"] = forcing_data["no3b"][month_counter]
    # month_data["po4b"] = forcing_data["po4b"][month_counter]
    # month_data["ponb"] = forcing_data["ponb"][month_counter]

    sclim = forcing_data[0]
    tclim = forcing_data[1]
    wclim = forcing_data[2]
    weddy1 = forcing_data[3]
    weddy2 = forcing_data[4]
    ism = forcing_data[5][:-1]
    wsu = -0.001 * forcing_data[6]
    wsv = -0.001 * forcing_data[7]
    swrad = -forcing_data[8] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
    wtsurf = -forcing_data[9] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
    qcorr = forcing_data[10]
    no3s = forcing_data[11]
    nh4s = forcing_data[12]
    po4s = forcing_data[13]
    sio4s = forcing_data[14]
    o2b = forcing_data[15]
    no3b = forcing_data[16]
    po4b = forcing_data[17]
    ponb = forcing_data[18]

    month_data = [sclim, tclim, wclim, weddy1, weddy2, ism, wsu, wsv, swrad, wtsurf, qcorr,
                  no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb]

    return month_data


# def forcing_manager(iter, physical, forcing, pom1d):
#     """
#     Description: Handles the forcing data reading and interpolation.
#                  Accomodates perpetual year monthly forcing time series.
#                  Data currently handled are time series of monthly averaged data.
                 
#     :return: forcing data
#     """
#     # INITIALISATION AND FIRST FORCING READING
#     if iter == 0:

#         # Month and day counters
#         forcing["counters"]["day_counter"] = 1
#         forcing["counters"]["month_counter"] = 0

#         # Timesteps to cover one day/month
#         forcing["counters"]["timesteps_per_day"] = physical["simulation"]["sec_per_day"] / physical["simulation"]["dt"]
#         forcing["counters"]["timesteps_per_month"] = 30 * forcing["counters"]["timesteps_per_day"]
        
#         # Month and day interpolators
#         forcing["counters"]["day_interpolator"] = -1
#         forcing["counters"]["month_interpolator"] = (forcing["counters"]["timesteps_per_month"] / 2.) - 1    # Climatological forcing data centered at day 15 of each month

#         # Reading for the first month
#         forcing["month1"] = write_forcing_data(forcing["counters"]["month_counter"], physical, pom1d)

#         # Update the day counter
#         forcing["counters"]["day_counter"] = forcing["counters"]["day_counter"] + 1

#         # Update the month counter
#         forcing["counters"]["month_counter"] = forcing["counters"]["month_counter"] + 1

#         # Reading for the second month
#         forcing["month2"] = write_forcing_data(forcing["counters"]["month_counter"], physical, pom1d)

#     # Update interpolation counters
#     forcing["counters"]["day_interpolator"] = forcing["counters"]["day_interpolator"] + 1
#     forcing["counters"]["day_ratio"] = forcing["counters"]["day_interpolator"] / forcing["counters"]["timesteps_per_day"]

#     if forcing["counters"]["day_ratio"] == 1:
#         x=1

#     forcing["counters"]["month_interpolator"] = forcing["counters"]["month_interpolator"] + 1
#     forcing["counters"]["month_ratio"] = forcing["counters"]["month_interpolator"] / forcing["counters"]["timesteps_per_month"]

#     # Interpolate wind stress
#     physical["stresses"]["wsu"] = forcing["month1"]["wsu"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["wsu"] - forcing["month1"]["wsu"])
#     physical["stresses"]["wsv"] = forcing["month1"]["wsv"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["wsv"] - forcing["month1"]["wsv"])

#     # Interpolate heat flux
#     if pom1d["general"]["idiagn"] == 0:
#         physical["temperature"]["surf_flux"] = forcing["month1"]["wtsurf"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["wtsurf"] - forcing["month1"]["wtsurf"])
#         physical["swrad"] = forcing["month1"]["swrad"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["swrad"] - forcing["month1"]["swrad"])
#     else:
#         physical["temperature"]["surf_flux"] = 0  # not needed for diagnostic mode (idiagn = 1), see 4.6.5 in manual
#         physical["swrad"] = forcing["month1"]["swrad"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["swrad"] - forcing["month1"]["swrad"])

#     # Interpolate temperature and salinity profiles
#     physical["temperature"]["ti"] = forcing["month1"]["tclim"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["tclim"] - forcing["month1"]["tclim"])
#     physical["salinity"]["si"] = forcing["month1"]["sclim"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["sclim"] - forcing["month1"]["sclim"])
#     physical["wgen"]  = forcing["month1"]["wclim"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["wclim"] - forcing["month1"]["wclim"])

#     if forcing["counters"]["month_ratio"] <= 0.5:
#         physical["weddy"] = forcing["month1"]["weddy1"]
#     else:
#         physical["weddy"] = forcing["month1"]["weddy2"]

#     if pom1d["general"]["idiagn"] == 0:
#         physical["temperature"]["surf"] = physical["temperature"]["ti"][0]
#         physical["salinity"]["surf"] = physical["salinity"]["si"][0]
#     elif pom1d["general"]["idiagn"] == 1:
#         physical["temperature"]["tf"] = physical["temperature"]["ti"]
#         physical["salinity"]["sf"] = physical["salinity"]["si"]

#     # Interpolate suspended inorganic matter
#     physical["ism"] = forcing["month1"]["ism"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["ism"] - forcing["month1"]["ism"])
    
#     # Interpolate surface and bottom nutrients
#     physical["nutrients"]["no3s"] = forcing["month1"]["no3s"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["no3s"] - forcing["month1"]["no3s"])
#     physical["nutrients"]["nh4s"] = forcing["month1"]["nh4s"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["nh4s"] - forcing["month1"]["nh4s"])
#     physical["nutrients"]["po4s"] = forcing["month1"]["po4s"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["po4s"] - forcing["month1"]["po4s"])
#     physical["nutrients"]["sio4s"] = forcing["month1"]["sio4s"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["sio4s"] - forcing["month1"]["sio4s"])

#     physical["nutrients"]["o2b"] = forcing["month1"]["o2b"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["o2b"] - forcing["month1"]["o2b"])
#     physical["nutrients"]["no3b"] = forcing["month1"]["no3b"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["no3b"] - forcing["month1"]["no3b"])
#     physical["nutrients"]["po4b"] = forcing["month1"]["po4b"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["po4b"] - forcing["month1"]["po4b"])
#     physical["nutrients"]["ponb_grad"] = forcing["month1"]["ponb"] + forcing["counters"]["month_ratio"] * (forcing["month2"]["ponb"] - forcing["month1"]["ponb"])

#     if forcing["counters"]["month_interpolator"] == forcing["counters"]["timesteps_per_month"]:

#         # Update the month counter
#         if os.path.exists('output_file.txt'):
#             output_file = open('output_file.txt','a')
#             output_file.write('month_counter = ')
#             output_file.write(str(forcing["counters"]["month_counter"]))
#             output_file.write('\n')
#             output_file.close()
#         else:
#             output_file = open('output_file.txt','w')
#             output_file.write('month_counter = ')
#             output_file.write(str(forcing["counters"]["month_counter"]))
#             output_file.write('\n')
#             output_file.close()

#         print('month_counter = ',forcing["counters"]["month_counter"])
#         forcing["counters"]["month_counter"] = forcing["counters"]["month_counter"] + 1
        
#         # Reset the interpolator
#         forcing["counters"]["month_interpolator"] = 0

#         # Shift the monthly data
#         forcing["month1"] = forcing["month2"]

#         # If 12 months have passed, restart the reading sequence
#         if forcing["counters"]["month_counter"] > 12:
#             forcing["counters"]["month_counter"] = 0
#             forcing["month1"] = write_forcing_data(forcing["counters"]["month_counter"], physical, pom1d)

#             forcing["counters"]["month_counter"] = forcing["counters"]["month_counter"] + 1

#         # Read the following month
#         forcing["month2"] = write_forcing_data(forcing["counters"]["month_counter"], physical, pom1d)

#     return physical, forcing



# def write_forcing_data(month_counter, physical, pom1d):
#     """
#     Definition: Writes forcing data for each month.
#     """
#     forcing_data = read_pom_inputs(physical, pom1d)
#     month_data = {}
    
#     # Climatology
#     month_data["sclim"] = forcing_data["sclim"][:,month_counter]
#     month_data["tclim"] = forcing_data["tclim"][:,month_counter]
#     month_data["wclim"] = forcing_data["wclim"][:,month_counter]

#     # Intermittant eddy w velocities
#     month_data["weddy1"] = forcing_data["weddy1"][:,month_counter]
#     month_data["weddy2"] = forcing_data["weddy2"][:,month_counter]

#     # Inorganic suspended matter
#     month_data["ism"] = forcing_data["ism"][:-1,month_counter]

#     month_data["wsu"] = -0.001 * forcing_data["wsu"][month_counter]     # N/m2-->m2/s2
#     month_data["wsv"] = -0.001 * forcing_data["wsv"][month_counter]     # N/m2-->m2/s2
    
#     # Heat Flux
#     month_data["swrad"]  = -forcing_data["swrad"][month_counter] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
#     month_data["wtsurf"] = -forcing_data["wtsurf"][month_counter] / pom1d["general"]["water_specific_heat_times_density"]   # W/m2-->deg.C*m/s
#     month_data["qcorr"]  = forcing_data["qcorr"][month_counter]

#     # Surface nutrients
#     month_data["no3s"]  = forcing_data["no3s"][month_counter]
#     month_data["nh4s"]  = forcing_data["nh4s"][month_counter]
#     month_data["po4s"]  = forcing_data["po4s"][month_counter]
#     month_data["sio4s"] = forcing_data["sio4s"][month_counter]

#     # Bottom nutrients
#     month_data["o2b"]  = forcing_data["o2b"][month_counter]
#     month_data["no3b"] = forcing_data["no3b"][month_counter]
#     month_data["po4b"] = forcing_data["po4b"][month_counter]
#     month_data["ponb"] = forcing_data["ponb"][month_counter]

#     return month_data
