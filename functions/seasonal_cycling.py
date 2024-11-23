import numpy as np

def calc_day_of_year(time_seconds):
    """ function that calculates the day of year """
    
    cycle = 366.0
    
    # calculate day of year
    if time_seconds == 0.0:
        day_of_year = 1.0
    else:
        sec_per_day = 86400.0
        time_day = time_seconds/sec_per_day
        day_of_year = np.floor(time_day) + 1
        
    # calculate year
    year = np.ceil(day_of_year/cycle)
    
    # correct day of year if year>1
    if year > 1:
        day_of_year -= (year - 1)*cycle

    return day_of_year

def calc_fraction_of_day(time_seconds):
    """ function that calculates the fration of the day """
    
    sec_per_day = 86400.0
    seconds_of_day = np.fmod(time_seconds,sec_per_day)
    fraction_of_day = seconds_of_day/sec_per_day

    return fraction_of_day


def get_wind(time,w_win,w_sum):
    """ function that calculates the seasonal wind values """

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    wind = (w_sum+w_win)/2 - ((w_sum-w_win)/2)*np.cos((day_of_year+(fraction_of_day - 0.5))*(np.pi/180))

    return wind


def get_salinity(time,s_win,s_sum):
    """ function that calculates the seasonal salinity values """

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    salinity = (s_sum+s_win)/2.0 - ((s_sum-s_win)/2.0)*np.cos((day_of_year+(fraction_of_day - 0.5))*(np.pi/180))

    return salinity


def get_mixed_layer_depth(time, mld_win, mld_sum):
    """ function that calculates the seasonal mixed layer depth values """

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    mld = (mld_sum+mld_win)/2.0 - ((mld_sum-mld_win)/2.0)*np.cos((day_of_year+(fraction_of_day - 0.5))*(np.pi/180))

    return mld


def get_sunlight(time,q_win,q_sum):
    """ function that calculates the sunlight """

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    latitude = 45.0
    light = (q_sum+q_win)/2.0 - (q_sum-q_win)/2.0*np.cos(day_of_year*(np.pi/180))
    cycle = 360
    declination = -0.406*np.cos(2.0*np.pi*int(day_of_year)/cycle)
    day_length = np.arccos(-np.tan(declination)*np.tan(latitude*(np.pi/180)))/np.pi*24.0
    day_time = fraction_of_day*24.0
    day_time = np.abs(day_time - 12.0)
    day_len = day_length/2.0
    if(day_time<day_len):
        day_time = day_time/day_len*np.pi
        wlight = light*np.cos(day_time) + light
    else:
        wlight = 0.0
    
    return wlight

def get_irrad(time,q_win,q_sum):

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    irrad = (q_sum + q_win)/2.0 - (q_sum - q_win)/2.0*np.cos((day_of_year+(fraction_of_day - 0.5))*(np.pi/180)) # - 0.5*np.cos(2*np.pi*fraction_of_day)

    return irrad

def get_temperature(time,t_win,t_sum):
    """ function that calculates the seasonal temperature """

    day_of_year = calc_day_of_year(time)
    fraction_of_day = calc_fraction_of_day(time)
    temperature = (t_sum + t_win)/2.0 - (t_sum - t_win)/2.0*np.cos((day_of_year+(fraction_of_day - 0.5))*(np.pi/180)) # - 0.5*np.cos(2*np.pi*fraction_of_day)

    return temperature