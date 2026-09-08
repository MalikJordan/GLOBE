import numpy as np
import sys
from numba import njit, types
from numba.typed import Dict, List
np.set_printoptions(precision=20)


def light_attenuation(base_element, light_attenuation_water, conc, tracer_map, tracers):
    """
    Definition:: Calculates light attenuation factor for photosynthesis
    Beer's Law attenuation coefficient
    """
    k_PAR = light_attenuation_water * np.ones(conc.shape[1],dtype=np.float64)

    for key in tracers:
        if tracers[key].type == "detritus":
            base_index = tracers[key].composition.index(base_element)
            k_PAR += tracers[key].light_attenuation * conc[tracer_map[key][base_index]]
        if tracers[key].type == "phytoplankton":
            # Extract light attenuation constant from parameter list
            lam = tracers[key].growth_ids.index("light_attenuation")
            light_coeff = tracers[key].growth_params[lam]

            # Light attenuation coefficient for phytoplankton is calculated using chl if available
            if "chl" in tracers[key].composition:
                chl_index = tracers[key].composition.index("chl")
                k_PAR += light_coeff * conc[tracer_map[key][chl_index]]
            else:
                base_index = tracers[key].composition.index(base_element)
                k_PAR += light_coeff * conc[tracer_map[key][base_index]]

    return k_PAR


@njit
def light_limitation(abbrev, growth_ids, growth_params, dz, irrad, k_PAR, Vm, temp_regulation_factor, nutrient_colimitation_factor, conc, tracer_map, composition):
    """
    k_PAR = Light Attenuation Coefficient
    surface_PAR = Photosynthetically Active Radiation (PAR) at watetr surface (z = 0)
    """
    light_limitation = growth_ids.index("light_limitation")
    # -------------------------------------------------------------------------------------------------
    # Monod
    # -------------------------------------------------------------------------------------------------
    if growth_params[light_limitation] == -1.:  # "monod"
        half_sat_light = growth_ids.index("half_sat_light")
        # Heinle & Slawig (2013)
        light_limitation_factor = irrad / (growth_params[half_sat_light] + irrad + 1E-20)

        # fake number
        irrad_at_depth = 1E-20 * np.ones(len(Vm))

    # -------------------------------------------------------------------------------------------------
    # Geider et al. (1997) / Jassby and Platt (1976)
    # -------------------------------------------------------------------------------------------------
    # elif growth_params[light_limitation] in [-2.,-3.]:  # ["geider","platt"]
    elif growth_params[light_limitation] == -2. or growth_params[light_limitation] == -3.:  # ["geider","platt"]
        light_location = growth_ids.index("light_location")
        initial_PI_slope = growth_ids.index("initial_PI_slope")

        # Calculate irradiance at depth
        if growth_params[light_location] == 1.:     # "top"
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * 86400                              # *86400 to convert from [1/s] to [1/d]
        elif growth_params[light_location] == 2.:   # "middle"  # Lazzari et al. (2012)
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * np.exp( -k_PAR * dz/2) * 86400     # *86400 to convert from [1/s] to [1/d]
        elif growth_params[light_location] == 3.:   # "integrated"  # Vichi et al. (2007)
            r = irrad / (k_PAR * dz) * (1. - np.exp(-k_PAR*dz))            
            irrad_at_depth = np.maximum(1E-20*np.ones_like(r), r*86400)                                        # *86400 to convert from [1/s] to [1/d]

        # Calculate Chl:C ratio using either ...   
        # if {"c","chl"}.issubset(composition): # Carbon and Chlorophyll concentrations (if preselt)
        if "c" in composition and "chl" in composition:
            carbon_index = composition.index("c")
            pc = conc[tracer_map[abbrev][carbon_index]]
            chl_index = composition.index("chl")
            pl = conc[tracer_map[abbrev][chl_index]]
            pl_pc = pl / pc     # Chl:C ratio (used in light limitation)
        else: # Geider et al. (1997) Dynamic Model
            theta_max = growth_ids.index("theta_max")
            if "theta_min" not in growth_ids: theta_min = 0.
            else:
                theta_min = growth_params[growth_ids.index("theta_min")]
            pl_pc = (growth_params[theta_max] - theta_min) / ( 1. + ( ( growth_params[theta_max] * growth_params[initial_PI_slope] * irrad_at_depth ) / 
                                                                                ( 2 * Vm * temp_regulation_factor * nutrient_colimitation_factor + 1.E-20 ) ) ) \
                        + theta_min

        # Calculate exponent for light limitation
        exp = pl_pc * ( growth_params[initial_PI_slope] / Vm ) * irrad_at_depth     # Stays like this for Platt
        if growth_params[light_limitation] == 2.:   # "geider"  # scale by temperature and nutrient limitation factors for Geider
            exp = exp / ( temp_regulation_factor * nutrient_colimitation_factor )

        light_limitation_factor = 1. - np.exp(-exp)

    # -------------------------------------------------------------------------------------------------
    # Smith (1936)
    # -------------------------------------------------------------------------------------------------
    elif growth_params[light_limitation] ==  -4.:   # "smith"
        light_location = growth_ids.index("light_location")
        initial_PI_slope = growth_ids.index("initial_PI_slope")
        
        # Calculate irradiance at depth
        if growth_params[light_location] == 1.:     # "top"
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad)
        elif growth_params[light_location] == 2.:   # "middle"  # Lazzari et al. (2012)
            # irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * np.exp( -k_PAR * dz/2)
            irrad_at_depth = np.maximum(1E-20*np.ones_like(irrad), irrad) * np.exp( -k_PAR * dz)
        elif growth_params[light_location] == 3.:   # "integrated"  # Vichi et al. (2007)
            r = irrad / (k_PAR * dz) * (1. - np.exp(-k_PAR*dz))            
            irrad_at_depth = np.maximum(1E-20*np.ones_like(r), r)  
        
        # Evans & Parslow (1985) formulation
        # num = Vm * growth_params[initial_PI_slope] *irrad_at_depth
        num = growth_params[initial_PI_slope] *irrad_at_depth
        den = np.sqrt((Vm**2) + ((growth_params[initial_PI_slope]*irrad_at_depth)**2))
        
        light_limitation_factor = num/(den + 1E-20)

    return irrad_at_depth, light_limitation_factor


@njit
def max_growth_rate(growth_ids, growth_params, temperature):
    """
    Defiition:: Calculates the temperature-dependent maximum phytoplankton growth rate, Vm
    """
    type = growth_ids.index("type")
    if growth_params[type] == 1.:    # "base_b"
        a = growth_ids.index("a")
        b = growth_ids.index("b")
        c = growth_ids.index("c")
        Vm = growth_params[a] * ( growth_params[b] ** ( growth_params[c] * temperature ) )

    elif growth_params[type] == 2.:  # "standard":
        base_growth_rate = growth_ids.index("a")
        eppley_coeff = growth_ids.index("eppley_coeff")
        Vm = growth_params[base_growth_rate] * np.exp(growth_params[eppley_coeff] * temperature)

    return Vm


@njit
def irradiance(eps_PAR, surface_PAR, depth, k_PAR):
    """
    Definition:: Calculates usable light for photosynthesis
    eps_PAR = fraction of photosynthetically available radiation
    0.217 = conversion from Einstein to Watts
    """

    # irradiance = surface_PAR * eps_PAR / 0.217

    irradiance = np.zeros(len(depth))
    irradiance[0] = surface_PAR * eps_PAR / 0.217
    if len(depth) > 1:
        for i in range(1,len(depth)):
            irradiance[i] = irradiance[i-1] * np.exp(-1. * k_PAR[i-1] * depth[i-1])

    return irradiance


@njit
def nutrient_limitation(nutrient, half_sat):
    """
    Definition:: Calculates the limitation factor a nutrient using the Michaelis-Menten formulation
    """
    nutrient_limitation_factor = nutrient / (half_sat + nutrient + 1.E-20)
    
    return nutrient_limitation_factor


@njit
def monod(nutrient, half_sat, exponent):

    limitation_factor = np.power(nutrient, exponent) / ( np.power(nutrient, exponent) + np.power(half_sat, exponent) + 1.E-20)

    return limitation_factor


@njit
def temperature_dependence(temperature, temp_reg_ids, temp_reg_params):
    """
    Definition:: Calculates temperature regulating factor
    """
    # Extract parameter indices
    func = temp_reg_ids.index("function")
    coeff = temp_reg_ids.index("coefficient")
    if "activation_energy" in temp_reg_ids: ae = temp_reg_ids.index("activation_energy")
    if "base_temp" in temp_reg_ids: bt = temp_reg_ids.index("base_temp")

    if temp_reg_params[func] == 1.:     # Arrhenius
        # Convert temperature from Celsius to Kelvin (+273.15)
        # Universal gas constant (R = 8.314 J mol-1 K-1)
        temp_regulating_factor = temp_reg_params[coeff] * np.exp(-temp_reg_params[ae] / (8.314 * (temperature+273.15)) )
    
    elif temp_reg_params[func] == 2.:   # Eppley
        temp_regulating_factor = np.exp(temp_reg_params[coeff] * temperature)
    
    elif temp_reg_params[func] == 3.:   # Q10
        # temp_regulating_factor = q10**((temperature-base_temp)/base_temp)
        temp_regulating_factor = temp_reg_params[coeff]**((temperature - temp_reg_params[bt]) / temp_reg_params[bt])

    return temp_regulating_factor


@njit
def concentration_ratio(conc, tracer_map):
    conc_ratio = np.ones_like(conc,dtype=np.float64)

    for tracer in tracer_map:
        if len(tracer_map[tracer]) > 1:
            base_conc = conc[tracer_map[tracer][0]]     # concentration of base element in composition
            
            for idx in range(1,len(tracer_map[tracer])):
                conc_ratio[tracer_map[tracer][idx]] = conc[tracer_map[tracer][idx]] / ( base_conc + 1.E-20 )
    
    return conc_ratio


def tracer_elements(base_element, reaction, tracers):
    """
    Definition:: Creates dictionary of tracer elements used for a particular reaction

    c/p   = consumed/produced tracer
    ic/ip = index of base element in consumed/produced tracer
    ec/ep = array of  elements in consumed/produced tracer affected by current reaction
    """
    if "consumed" in reaction and reaction["consumed"] != None:    rc = reaction["consumed"]
    else:   rc = []
    if "produced" in reaction and reaction["produced"] != None:    rp = reaction["produced"]
    else:   rp = []

    
    # Initialize numba typed.Dicts
    ic = Dict.empty(key_type=types.unicode_type, value_type=types.int64[:])
    ec = Dict.empty(key_type=types.unicode_type, value_type=types.int64[:])
    ip = Dict.empty(key_type=types.unicode_type, value_type=types.int64[:])
    ep = Dict.empty(key_type=types.unicode_type, value_type=types.int64[:])

    if rc:
        for key in rc:
            if tracers[key].type == "inorganic":
                ic[key] = np.array([0],dtype=np.int64)
                ec[key] = np.array([1],dtype=np.int64)
            else:
                ic[key] = np.array([tracers[key].composition.index(base_element)],dtype=np.int64)
                ec[key] = np.zeros(len(tracers[key].composition),dtype=np.int64)
                for element in tracers[key].composition:
                    i = tracers[key].composition.index(element)
                    # if element in reaction["consumed"][key]:  ec[key][i] = np.array([1],dtype=np.int64)
                    if element in reaction["consumed"][key]:  ec[key][i] = np.int64(1)
    
    if rp:
        for key in rp:
            if tracers[key].type == "inorganic":
                ip[key] = np.array([0],dtype=np.int64)
                ep[key] = np.array([1],dtype=np.int64)
            else:
                ip[key] = np.array([tracers[key].composition.index(base_element)],dtype=np.int64)
                ep[key] = np.zeros(len(tracers[key].composition),dtype=np.int64)
                for element in tracers[key].composition:
                    i = tracers[key].composition.index(element)
                    # if element in reaction["produced"][key]:  ep[key][i] = np.array([1],dtype=np.int64)
                    if element in reaction["produced"][key]:  ep[key][i] = np.int64(1)

    c = List(rc.keys()) if rc else List.empty_list(types.unicode_type)
    p = List(rp.keys()) if rp else List.empty_list(types.unicode_type)

    return c, p, ec, ep, ic, ip


@njit
def switch(parameter):

    x = len(parameter)
    if x > 1:
        x = np.shape(parameter)[0]
        switch = np.zeros(x)
        for i in range(0,x):
            if parameter[i] > 0.0:
                switch[i] = 1.0
    else:
        if parameter > 0.:  switch = np.array([1.],dtype=np.float64)
        else:   switch = np.array([0.],dtype=np.float64)

    return switch


@njit
def string_to_float(s):
    if len(s) == 0:
        return 0.0

    value = 0.0
    decimal_place = 0.0
    exponent = 0
    exponent_sign = 1

    sign = 1.0
    i = 0

    # leading sign
    if s[0] == '-':
        sign = -1.0
        i += 1
    elif s[0] == '+':
        i += 1

    if i == len(s):
        return 0.0

    has_digit = False
    has_decimal = False
    has_exponent = False

    # mantissa
    while i < len(s):
        char = s[i]

        if char >= '0' and char <= '9':
            has_digit = True
            digit = ord(char) - ord('0')

            if not has_decimal:
                value = value * 10.0 + digit
            else:
                decimal_place += 1.0
                value += digit / (10.0 ** decimal_place)

        elif char == '.':
            if has_decimal or has_exponent:
                return 0.0
            has_decimal = True

        elif char == 'e' or char == 'E':
            has_exponent = True
            i += 1
            break

        else:
            return 0.0

        i += 1

    if not has_digit:
        return 0.0

    # exponent
    if has_exponent:
        if i == len(s):
            return 0.0

        if s[i] == '-':
            exponent_sign = -1
            i += 1
        elif s[i] == '+':
            i += 1

        if i == len(s):
            return 0.0

        exp_value = 0
        exp_digit = False

        while i < len(s):
            char = s[i]

            if char >= '0' and char <= '9':
                exp_digit = True
                exp_value = exp_value * 10 + (ord(char) - ord('0'))
            else:
                return 0.0

            i += 1

        if not exp_digit:
            return 0.0

        exponent = exponent_sign * exp_value

    return sign * value * (10.0 ** exponent)


@njit
def calculate_acidity(temperature, salinity, density, wind, d, schmidt_ratio, co2, ta):
    """ This function calculates the solubility and acidity constants (k values).
    """
    # Compute Chemical enhancement the Temperature dependent gas transfer
    bt = 2.5*(0.5246 + 1.6256e-2*temperature + 4.9946e-4*(temperature**2))

    # Calculate wind dependency + Chemical enhancement including conversion cm/hr => m/s (0.01 = cm to m, 24 = days to hours)
    ken = (bt + d*(wind**2))* np.sqrt(schmidt_ratio)* 0.01 * 24.

    # K0, solubility of co2 in the water (K Henry) from Weiss 1974; K0 = [co2]/pco2 [mol kg-1 atm-1]
    # 273.15 = Celsius to Kelvin
    tk = (temperature + 273.15)
    tk_100 = tk/100.0
    k0 = np.exp((93.4517/tk_100) - 60.2409 + (23.3585*np.log(tk_100)) + salinity*(0.023517 - (0.023656*tk_100) + 0.0047036*(tk_100**2)))
    
    # calculate all constants needed to convert between various measured carbon species
    tk_inv = 1.0/tk
    tk_log = np.log(tk)
    salt_sqrt = np.sqrt(salinity)
    
    # chlorinity
    scl = salinity/1.80655
    
    # ionic strength 
    ionic_strength = 19.924*salinity/(1000.0 - 1.005*salinity)
    ionic_strength_sqrt = np.sqrt(ionic_strength)
    
    # Calculate concentrations for borate, sulfate, and fluoride as a function of chlorinity
    # Uppstrom (1974)
    bt = 0.000232*scl/10.811
    
    # Morris & Riley (1966)
    st = 0.14*scl/96.062
    
    # Riley (1965)
    ft = 0.000067*scl/18.9984
    
    # change units from mmol m^-3 to mol/kg
    pt = 0.0
    sit = 0.0
    
    # change DIC (o3h) units from mg C m^-3 to mol kg^-1
    # 0.001 = g/mg
    # 12 = stoichiometric coefficient between Carbon and O2 [C:O2 = 1:12]
    ldic = co2 * (1/12) * 0.001 / density
    
    # convert TA from inits of mmol eq m^-3 to mol eq kg^-1
    # 0.001 = mol/mmol
    alk = ta / density * 0.001
    
    # calculate Acidity constants
    # constants according to Mehrbach et al. (1973) as refitted by Dickson and Millero (1987) 
    # ph scale: seawater  (Millero, 1995, p.664)
    # Standard OCMIP computation. Natural seawater
    # calculate carbonate equilibrium I, k1, [H][hco3]/[H2co3]
    k1 = 10.0**(-1.0*(3670.7*tk_inv - 62.008 + 9.7944*tk_log - 0.0118*salinity + 0.000116*(salinity**2)))

    # calculate carbonate equilibrium II, k2, [H][co3]/[hco3]
    k2 = 10.0**(-1.0*(1394.7*tk_inv + 4.777 - 0.0184*salinity + 0.000118*(salinity**2)))

    # calculate k1p, [H][H2PO4]/[H3PO4] (ph scale: total)
    lnK = -4576.752*tk_inv + 115.525 - 18.453* tk_log + (-106.736*tk_inv + 0.69171)*salt_sqrt + (-0.65643*tk_inv - 0.01844)*salinity
    k1p = np.exp(lnK)

    # calculate k2p, [H][HPO4]/[H2PO4] (ph scale: total)
    lnK = -8814.715*tk_inv + 172.0883 - 27.927* tk_log +(-160.340*tk_inv + 1.3566)*salt_sqrt + (0.37335*tk_inv - 0.05778)*salinity
    k2p = np.exp(lnK)
    
    # calculate k3p, [H][PO4]/[HPO4] (ph scale: total)
    lnK = -3070.75*tk_inv - 18.126 + (17.27039*tk_inv + 2.81197)*salt_sqrt + (-44.99486*tk_inv - 0.09984)*salinity
    k3p = np.exp(lnK)
    
    # calculate ksi, [H][SiO(OH)3]/[Si(OH)4]
    ksi = -8904.2*tk_inv + 117.385 - 19.334*tk_log + (-458.79*tk_inv + 3.5913)*ionic_strength_sqrt + (188.74*tk_inv - 1.5998)*ionic_strength + (-12.1652*tk_inv + 0.07871)*(ionic_strength**2) + np.log(1.0 - (0.001005*salinity))
    ksi = np.exp(lnK)
    
    # calculate kw, [H][OH] (ph scale: SWS)
    intercept = 148.9802
    lnK = intercept - 13847.26*tk_inv - 23.6521*tk_log + (118.67*tk_inv - 5.977 + 1.0495*tk_log)*salt_sqrt - 0.01615*salinity
    kw = np.exp(lnK)
    
    # calculate ks, [H][SO4]/[HSO4]  (ph scale: "free")
    lnK = -4276.1*tk_inv + 141.328 - 23.093*tk_log + (-13856.0*tk_inv + 324.57 - 47.986*tk_log)*ionic_strength_sqrt + (35474.0*tk_inv - 771.54 + 114.723*tk_log)*ionic_strength - 2698.0*tk_inv*(ionic_strength**1.5) + 1776.0*tk_inv*(ionic_strength**2) + np.log(1.0 - 0.001005*salinity)
    ks = np.exp(lnK)
    
    # calculate kf, [H][F]/[HF] (ph scale: "free")
    lnK = 1590.2*tk_inv - 12.641 + 1.525*ionic_strength_sqrt + np.log(1.0 - 0.001005*salinity)
    kf = np.exp(lnK)
    
    # calculate kb, [H][BO2]/[HBO2] (ph scale: total)
    lnK = (-8966.90 - 2890.53*salt_sqrt - 77.942*salinity + 1.728*(salinity**1.5) - 0.0996*(salinity**2))*tk_inv + (148.0248 + 137.1942*salt_sqrt + 1.62142*salinity) + (-24.4344 - 25.085*salt_sqrt - 0.2474*salinity)*tk_log + 0.053105*salt_sqrt*tk
    kb = np.exp(lnK)

    return k0, k1, k1p, k2, k2p, k3p, ksi, kw, ks, kf, kb, ken, bt, st, ft, pt, sit, ldic, alk


# @njit
def calculate_Hplus(pH, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk):
    """ This function expresses total alkalinity (TA) as a function of DIC, 
    hSWS (H+ on sea water scale) and constants. It also calculates the 
    derivative of this function with respect to hSWS.
    
    This function was obtained from ModuleCO2System.F90
    
    Function inputs:
        x:           H+ on sea water scale
    
    Function outputs:
        f_TA:        calculated value for TA
        df_TA_dhSWS: derivative of f_TA with respect to hSWS
    """
    
    t1 = 1.0
    t2 = 2.0
    t3 = 3.0
    
    # derive H+ and other constants
    pH_2 = pH*pH
    pH_3 = pH_2*pH
    k12 = k1*k2
    k12p = k1p*k2p
    k123p = k12p*k3p
    c = t1 + st/ks + ft/kf
    a = pH_3 + k1p*pH_2 + k12p*pH + k123p
    a2 = a*a
    da = t3*pH_2 + t2*k1p*pH + k12p
    b = pH_2 + k1*pH + k12
    b2 = b*b
    db = t2*pH + k1
    
    fn = k1*pH*ldic/b + t2*ldic*k12/b + bt/(t1 + pH/kb)
    df = ((k1*ldic*b) - k1*pH*ldic*db)/b2 - t2*ldic*k12*db/b2 - bt/kb/(t1 + pH/kb)**t2
    
    fn += (kw/pH + pt*k12p*pH/a + t2*pt*k123p/a + sit/(t1 + pH/ksi) - pH/c - 
          st/(t1 + ks/(pH/c)) - ft/(t1 + kf/(pH/c)) - pt*pH_3/a - alk)
    
    df += (-kw/pH_2 + (pt*k12p*(a - pH*da))/a2 - t2*pt*k123p*da/a2 - 
          sit/ksi/(t1 + pH/ksi)**t2 - t1/c - st*(t1 + ks/(pH/c))**(-t2)*(ks*c/pH_2) - 
          ft*(t1 + kf/(pH/c))**(-t2)*(kf*c/pH_2) - pt*pH_2*(t3*a - pH*da)/a2)
    
    return fn, df


# @njit
def find_roots_of_f_TA(x1, x2, xacc, maxit, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk):
    """ This function finds the roots of the total alkalinity function
    
    This function was obtained from ModuleCO2System.F90
    
    Function inputs:
        x1
        x2
        xacc:           Accuracy of the iterative scheme for OCMIP
        maxit:          Maximum number of iterations for OCMIP
    
    Function outputs:
        error
    """
    
    error = 0
    fl, df = calculate_Hplus(x1, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
    fh, df = calculate_Hplus(x2, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
    if fl==0:
        drtsafe2 = x1
        error = 1
    elif fh ==0:
        drtsafe2 = x2
        error = 1
    elif fl<0:
        xl = x1
        xh = x2
    else:
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap

    drtsafe2 = 0.5*(x1 + x2)
    dxold = np.abs(x2 - x1)
    dx = dxold
    f, df = calculate_Hplus(drtsafe2, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
    
    j = 0
    ready = False
    
    while True:
        j+=1
        if ((drtsafe2 - xh)*df - f)*((drtsafe2 - xl)*df - f) >= 0 or np.abs(2.0*f) > np.abs(dxold*df):
            dxold = dx
            dx = 0.5*(xh - xl)
            drtsafe2 = xl + dx
            ready = (xl == drtsafe2)
        else:
            dxold = dx
            dx = f/df
            temp = drtsafe2
            drtsafe2 = drtsafe2 - dx
            ready = (temp == drtsafe2)
        ready = np.abs(dx)<xacc
        if not ready:
            f, df = calculate_Hplus(drtsafe2, k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
            if f<0:
                xl = drtsafe2
                fl = f
            else:
                xh = drtsafe2
                fh = f
        if ready and j<maxit:
            break
        if not ready and np.isnan(drtsafe2):
            drtsafe2 = 10**(8.12)     # pH = log10(drtsafe2) set to initial value provided in co2_flux_parameters
            break
    if j>maxit:
        error = 2
    
    return drtsafe2, error
