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
    sclim = forcing_data[0]
    tclim = forcing_data[1]
    wclim = forcing_data[2]

    # Intermittant eddy w velocities
    weddy1 = forcing_data[3]
    weddy2 = forcing_data[4]

    # Inorganic suspended matter
    ism = forcing_data[5][:-1]

    # Wind stress
    wsu = -0.001 * forcing_data[6]
    wsv = -0.001 * forcing_data[7]

    # Heat Flux
    swrad = -forcing_data[8] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
    wtsurf = -forcing_data[9] / pom1d["general"]["water_specific_heat_times_density"]    # W/m2-->deg.C*m/s
    qcorr = forcing_data[10]

    # Surface nutrients
    no3s = forcing_data[11]
    nh4s = forcing_data[12]
    po4s = forcing_data[13]
    sio4s = forcing_data[14]

    # Bottom nutrients
    o2b = forcing_data[15]
    no3b = forcing_data[16]
    po4b = forcing_data[17]
    ponb = forcing_data[18]

    month_data = [sclim, tclim, wclim, weddy1, weddy2, ism, wsu, wsv, swrad, wtsurf, qcorr,
                  no3s, nh4s, po4s, sio4s, o2b, no3b, po4b, ponb]

    return month_data
