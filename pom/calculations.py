import numpy as np
from numba import njit
from pom.check_diffusion import a, c, vh, vhp
np.set_printoptions(precision=20)

# order of inputs --> (case, dt2, num_layers, column_depth, z, dz, zz, dzz, l, mld, umol, nbc, ntp, swrad, km, kh, kq, ke, kef, keb, kel, kelf, kelb, density, temp, sal, u, uf, ub, bsu, wsu, v, vf, vb, bsv, wsv)
# @njit
# def density_profile(num_layers, column_depth, dzz, temp, sal):
#     """
#     Description: Calculates vertical density profile
                 
#                  dzz = staggered vertical spacing
#     """
    
#     gravity = 9.806
#     vertical_density_profile = np.zeros(num_layers)
#     pressure = -gravity * 1.025 * dzz[:-1] * column_depth * 0.01

#     cr = 1449.1 + 0.0821*pressure + (4.55*temp[:-1]) - (0.045*np.power(temp[:-1],2)) + (1.34*(sal[:-1] - 35.))
#     cr = pressure/np.power(cr,2)

#     density = 999.842594 + (6.793952E-02 * temp[:-1]) - (9.095290E-03 * np.power(temp[:-1],2)) \
#             + (1.001685E-04 * np.power(temp[:-1],3)) - (1.120083E-06 * np.power(temp[:-1],4)) \
#             + (6.536332E-09 * np.power(temp[:-1],5))
#     density = density + (0.824493 - (4.0899E-03 * temp[:-1]) + (7.6438E-05 * np.power(temp[:-1],2)) \
#             - (8.2467E-07 * np.power(temp[:-1],3)) + (5.3875E-09 * np.power(temp[:-1],4))) * sal[:-1] \
#             + (-5.72466E-03 + (1.0227E-04 * temp[:-1]) - (1.6546E-06 * np.power(temp[:-1],2))) * (np.power(np.abs(sal[:-1]),1.5)) \
#             + 4.8314E-4*np.power(sal[:-1],2)
#     density = density + 1.E05*cr*(1. - (2*cr))

#     vertical_density_profile[:-1] = (density - 1000.) * 1.E-03
#     vertical_density_profile[-1] = vertical_density_profile[-2]

#     return vertical_density_profile


def density_profile(configuration, num_layers, column_depth, dzz, temp, sal):
    """
    Description: Calculates vertical density profile
                 
                 dzz = staggered vertical spacing
    """
    
    gravity = 9.806
    vertical_density_profile = np.zeros(num_layers)
    if configuration == "0d":
        pressure = -gravity * 1.025 * column_depth * 0.01
        
        cr = 1449.1 + 0.0821*pressure + (4.55*temp) - (0.045*np.power(temp,2)) + (1.34*(sal - 35.))
        cr = pressure/np.power(cr,2)

        density = 999.842594 + (6.793952E-02 * temp) - (9.095290E-03 * np.power(temp,2)) \
                + (1.001685E-04 * np.power(temp,3)) - (1.120083E-06 * np.power(temp,4)) \
                + (6.536332E-09 * np.power(temp,5))
        density = density + (0.824493 - (4.0899E-03 * temp) + (7.6438E-05 * np.power(temp,2)) \
                - (8.2467E-07 * np.power(temp,3)) + (5.3875E-09 * np.power(temp,4))) * sal \
                + (-5.72466E-03 + (1.0227E-04 * temp) - (1.6546E-06 * np.power(temp,2))) * (np.power(np.abs(sal),1.5)) \
                + 4.8314E-4*np.power(sal,2)
        density = density + 1.E05*cr*(1. - (2*cr))

        vertical_density_profile = (density - 1000.) * 1.E-03
        
    else:
        pressure = -gravity * 1.025 * dzz[:-1] * column_depth * 0.01

        cr = 1449.1 + 0.0821*pressure + (4.55*temp[:-1]) - (0.045*np.power(temp[:-1],2)) + (1.34*(sal[:-1] - 35.))
        cr = pressure/np.power(cr,2)

        density = 999.842594 + (6.793952E-02 * temp[:-1]) - (9.095290E-03 * np.power(temp[:-1],2)) \
                + (1.001685E-04 * np.power(temp[:-1],3)) - (1.120083E-06 * np.power(temp[:-1],4)) \
                + (6.536332E-09 * np.power(temp[:-1],5))
        density = density + (0.824493 - (4.0899E-03 * temp[:-1]) + (7.6438E-05 * np.power(temp[:-1],2)) \
                - (8.2467E-07 * np.power(temp[:-1],3)) + (5.3875E-09 * np.power(temp[:-1],4))) * sal[:-1] \
                + (-5.72466E-03 + (1.0227E-04 * temp[:-1]) - (1.6546E-06 * np.power(temp[:-1],2))) * (np.power(np.abs(sal[:-1]),1.5)) \
                + 4.8314E-4*np.power(sal[:-1],2)
        # density = density + 1.E05*cr*(1. - (2*cr))

        vertical_density_profile[:-1] = (density - 1000.) * 1.E-03
        vertical_density_profile[-1] = vertical_density_profile[-2]

    return vertical_density_profile


@njit
def kinetic_energy_profile(dt2, num_layers, column_depth, z, dz, dzz, l, umol, km, kh, kq, ke, kef, keb, kel, kelf, kelb, density, u, bsu, wsu, v, bsv, wsv):
    """
    Description: Solves for turbulent closure
                 
                 dt2 = twice the timestep
                 z = vertical grid
                 dz = vertical spacing
                 dzz = staggered vertical spacing
                 l = length scale
                 umol = background diffusion coefficient
                 
                 km = diffusion profile for momentum
                 kh = diffusion profile for tracers
                 kq = diffusion profile for kinetic energy

                 ke = kinetic energy at current step
                 kef = "                " forward step
                 keb = "                " backward step

                 kel = kinetic energy times length scale at current step
                 kelf = "                                   " forward step
                 kelb = "                                   " backward step

                 u = zonal velocity
                 wsu = zonal wind stress
                 bsu = zonal bottom stress
                 
                 v = meridional velocity
                 wsv = meridional wind stress
                 bsv = meridional bottom stress

    :return: data arrays for kinetic energy & kinetic energy times length profiles at foward step (kef,kelf), length scale, and diffusion profiles (km, kh, kq)
    """
    
    # INITIALIZE VARIABLES
    A = np.zeros(num_layers)
    C = np.zeros(num_layers)
    VH = np.zeros(num_layers)
    VHP = np.zeros(num_layers)

    KN = np.zeros(num_layers)
    GH = np.zeros(num_layers)
    SH = np.zeros(num_layers)
    SM = np.zeros(num_layers)

    DTEF = np.zeros(num_layers)
    BPROD = np.zeros(num_layers)
    PROD = np.zeros(num_layers)
    SPROD = np.zeros(num_layers)
    pressure = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOYGR = np.zeros(num_layers); CC = np.zeros(num_layers)
    TEMP1 = np.zeros(num_layers); TEMP2 = np.zeros(num_layers); TEMP3 = np.zeros(num_layers)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   DATA STATEMENTS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    A1 = 0.92
    B1 = 16.6
    A2 = 0.74
    B2 = 10.1
    C1 = 0.08
    E1 = 1.8
    E2 = 1.33
    E3 = 1.0
    von_karman_constant = 0.40  # KAPPA
    SQ = 0.2
    CIWC = 1.0
    gravity = 9.806
    SMALL = 1.E-08
    
    A[1:-1] = -dt2 * (kq[2:] + kq[1:-1] + 2 * umol) * \
                0.5 / (dzz[:-2] * dz[1:-1] * column_depth * column_depth)
    C[1:-1] = -dt2 * (kq[:-2] + kq[1:-1] + 2 * umol) * \
                0.5 / (dzz[:-2] * dz[:-2] * column_depth * column_depth)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR THE EQUATION
    #   DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    CONST1 = 16.6 ** 0.6666667 * CIWC

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    VH[0] = 0.0
    VHP[0] = np.sqrt(wsu ** 2 + wsv ** 2) * CONST1
    kef[num_layers - 1] = 0.5 * np.sqrt((bsu + bsu) ** 2 + (bsv + bsv) ** 2) * CONST1

    keb[1:-1] = np.abs(keb[1:-1])
    kelb[1:-1] = np.abs(kelb[1:-1])
    BOYGR[1:-1] = gravity * (density[:-2] - density[1:-1]) / (dzz[:-2] * column_depth)
    DTEF[1:-1] = keb[1:-1] * np.sqrt(keb[1:-1]) / (B1 * kelb[1:-1] + SMALL)
    SPROD[1:-1] = .25 * km[1:-1] * ((u[1:-1] + u[1:-1] - u[0:-2] - u[0:-2]) ** 2 + (v[1:-1] + v[1:-1] - v[0:-2] - v[0:-2]) ** 2) / (dzz[0:-2] * column_depth) ** 2 * CIWC ** 2
    BPROD[1:-1] = kh[1:-1] * BOYGR[1:-1]
    PROD[1:-1] = SPROD[1:-1] + BPROD[1:-1]
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for i in range(1, num_layers - 1):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (2. * dt2 * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (-2. * dt2 * PROD[i] + C[i] * VHP[i - 1] - keb[i]) * VHP[i]
    
    x=1
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for i in range(2, num_layers+1):  # 104
        k = num_layers - i
        kef[k] = VH[k] * kef[k + 1] + VHP[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    kelf[num_layers - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, num_layers - 1):
        DTEF[i] = DTEF[i] * (1. + E2 * ((1. / np.abs(z[i] - z[0]) + 1. / np.abs(z[i] - z[num_layers-1])) * l[i] / (column_depth * von_karman_constant)) ** 2)
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (dt2 * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (dt2 * (- (SPROD[i] + E3 * BPROD[i]) * l[i] * E1) + C[i] * VHP[i - 1] - kelb[i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(2, num_layers+1):
        k = num_layers - i
        kelf[k] = VH[k] * kelf[k + 1] + VHP[k]

    for i in range(1, num_layers - 1):
        if kef[i] > SMALL or kelf[i] > SMALL:
            continue
        else:
            kef[i] = SMALL
            kelf[i] = SMALL

    kef[:-1] = np.abs(kef[:-1])
    kelf[:-1] = np.abs(kelf[:-1])
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR KM AND KH
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    COEF1 = A2 * (1. - 6. * A1 / B1)
    COEF2 = 3. * A2 * B2 + 18. * A1 * A2
    COEF3 = A1 * (1. - 3. * C1 - 6. * A1 / B1)
    COEF4 = 18. * A1 * A1 + 9. * A1 * A2
    COEF5 = 9. * A1 * A2

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NOTE THAT SM AND SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    l[0] = 0.
    l[num_layers - 1] = 0.
    GH[0] = 0.
    GH[num_layers - 1] = 0.

    l[1:-1] = kelf[1:-1] / kef[1:-1]
    GH[1:-1] = np.power(l[1:-1],2) / kef[1:-1] * BOYGR[1:-1]
    for i in range(0, num_layers):
        GH[i] = min(GH[i], .028)
    SH = COEF1 / (1. - COEF2 * GH)
    SM = COEF3 + SH * COEF4 * GH
    SM = SM / (1. - COEF5 * GH)

    KN = l * np.sqrt(np.abs(ke))
    kq = 0.5 * (KN * 0.41 * SM + kq)
    km = 0.5 * (KN * SM + km)
    kh = 0.5 * (KN * SH + kh)

    return kef, kelf, l, km, kh, kq


@njit
def meridional_velocity_profile(dt2, num_layers, column_depth, dz, dzz, umol, km, vf, bsv, wsv):
    """ 
    Description: Calculates meridional (V) velocity profile
                 Solves for the equation:    dti2 * (KM * V')' - V = -VB
                 
                 dt2 = twice the timestep
                 dz = vertical spacing
                 dzz = staggered vertical spacing
                 km = diffusion of momentum
                 vf = forward step for meridional velocity
                 bsv = meridional bottom stress
                 wsv = meridional wind stress
                 umol = background diffusion coefficient
    
    :return: data array for meridional velocity profile
    """
    A = np.zeros(num_layers)
    C = np.zeros(num_layers)
    VH = np.zeros(num_layers)
    VHP = np.zeros(num_layers)
  
    A[:-2] = -dt2 * (km[1:-1] + umol) / (dz[:-2] * dzz[:-2] * column_depth * column_depth)
    C[1:-1] = -dt2 * (km[1:-1] + umol) / (dz[1:-1] * dzz[:-2] * column_depth * column_depth)

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-dt2 * wsv / (-dz[0] * column_depth) - vf[0]) / (A[0] - 1.)

    for i in range(1, num_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - vf[i]) * VHP[i]

    CBC = 0.0
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TO RESTORE BOTTOM B.L. DELETE NEXT LINE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    vf[num_layers - 2] = (C[num_layers - 1] * VHP[num_layers - 3] - vf[num_layers - 2]) / (CBC * dt2 / (-dz[num_layers - 2] * column_depth) - 1. - (VH[num_layers - 3] - 1.) * C[num_layers - 2])

    for i in range(1, num_layers - 1):
        k = num_layers - 1 - i
        vf[k - 1] = VH[k - 1] * vf[k] + VHP[k - 1]

    bsv = -CBC * vf[num_layers - 2]

    return vf, bsv


@njit
def mixed_layer_depth(num_layers, zz, mld, temp):
    """
    Description: Calculates meridional (V) velocity profile
                 Solves for the equation:    dti2 * (KM * V')' - V = -VB
                 
                 mld = mixed layer depth
                 zz = staggered vertical grid

    :return: data array for mixed layer depth
    """
    small = 1.e-06
    
    for i in range(0,num_layers-1):

        if temp[0] > temp[i]+0.2:
            break

        mld[i] = zz[i] - (temp[i] + 0.2 - temp[0]) * (zz[i] - zz[i+1]) / (temp[i] - temp[i+1] + small)
        
    return mld


@njit
def temperature_and_salinity_profiles(case, dt2, num_layers, column_depth, z, dz, dzz, umol, nbc, ntp, swrad, kh, forward, surface_value, surface_flux, bottom_flux):
    """ 
    Description: Solves for the conservative (temperature & salinity) and non-conservative (BFM state variables) scalars
                 Handles the surface and bottom boundary conditions
    NOTE: Conservative scalars are only calculated when the system is run in prognostic mode
    
    :return: data arrays for the 'property' (temperature, salinity, or BFM state variable)
    """
    
    # FLAG FOR BOUNDARY CONDITION DEFINITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   NBC=1: SURF. B.C. IS WFSURF+SWRAD. NO RADIATIVE PENETRATION.
    #   NBC=2; SURF. B.C. IS WFSURF. SWRAD PENETRATES WATER COLUMN.
    #   NBC=3; SURF. B.C. IS TSURF. NO SWRAD RADIATIVE PENETRATION.
    #   NBC=4; SURF. B.C. IS TSURF. SWRAD PENETRATES WATER COLUMN.
    #
    #   NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEGATIVE VALUES WHEN FLUX IS "IN" THE WATER COLUMN
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # FLAG FOR JERLOV WATER TYPE CHOICE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   JERLOV WATER TYPE CHOICE IS RELEVANT ONLY WHEN NBC = 2 OR NBC = 4.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    A = np.zeros(num_layers)
    C = np.zeros(num_layers)
    VH = np.zeros(num_layers)
    VHP = np.zeros(num_layers)

    # SW PROFILE
    vertical_radiation_profile = np.zeros(num_layers)

    # IRRADIANCE PARAMETERS AFTER PAULSON & SIMPSON JPO 1977, 952-956
    RP = [0.58, 0.62, 0.67, 0.77, 0.78]
    AD1 = [0.35, 0.60, 1.00, 1.50, 1.40]
    AD2 = [23.00, 20.00, 17.00, 14.00, 7.90]

    # JERLOV WATER TYPES
    # NTP         = 1           2            3           4          5
    # JERLOV TYPE = I           IA           IB          II         III
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   START COMPUTATION OF VERTICAL PROFILE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    A[:-2] = -dt2 * (kh[1:-1] + umol) / (dz[:-2] * dzz[:-2] * column_depth * column_depth)
    C[1:-1] = -dt2 * (kh[1:-1] + umol) / (dz[1:-1] * dzz[:-2] * column_depth * column_depth)

    vertical_radiation_profile[:] = 0.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.
    
    if nbc == 1:
        VH[0] = A[0] / (A[0] - 1.)
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            VHP[0] = -dt2 * (surface_flux + 0.) / (-dz[0] * column_depth) - forward[0]
        else:
            VHP[0] = -dt2 * (surface_flux + swrad) / (-dz[0] * column_depth) - forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 2:
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            vertical_radiation_profile[:] = 0. * (RP[ntp] * np.exp(z[:] * column_depth / AD1[ntp]) + (1. - RP[ntp] * np.exp(z[:] * column_depth / AD2[ntp])))  # ***
        else:
            vertical_radiation_profile[:] = swrad * (RP[ntp] * np.exp(z[:] * column_depth / AD1[ntp]) + (1. - RP[ntp] * np.exp(z[:] * column_depth / AD2[ntp])))  # ***
        vertical_radiation_profile[num_layers - 1] = 0.

        VH[0] = A[0] / (A[0] - 1.)
        VHP[0] = dt2 * (surface_flux + vertical_radiation_profile[0] - vertical_radiation_profile[1]) / (dz[0] * column_depth) - forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 3:
        VH[0] = 0.
        VHP[0] = surface_value

    elif nbc == 4:
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            vertical_radiation_profile[:] = 0. * (RP[ntp] * np.exp(z[:] * column_depth / AD1[ntp]) + (1. - RP[ntp] * np.exp(z[:] * column_depth / AD2[ntp])))  # ***
        else:
            vertical_radiation_profile[:] = swrad * (RP[ntp] * np.exp(z[:] * column_depth / AD1[ntp]) + (1. - RP[ntp] * np.exp(z[:] * column_depth / AD2[ntp])))  # ***
        vertical_radiation_profile[num_layers - 1] = 0.

        VH[0] = 0.
        VHP[0] = surface_value

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, num_layers - 2):
        VHP[i] = 1 / (A[i] + C[i] * (1 - VH[i - 1]) - 1)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - forward[i] + dt2 * (vertical_radiation_profile[i] - vertical_radiation_profile[i + 1]) / (column_depth * dz[i])) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    forward[num_layers - 2] = (C[num_layers - 2] * VHP[num_layers - 3] - forward[num_layers - 2] + (bottom_flux * dt2 / (dz[num_layers - 2] * column_depth))) / (C[num_layers - 2] * (1 - VH[num_layers - 3]) - 1)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(3, num_layers+1):
        k = (num_layers) - i
        forward[k] = VH[k] * forward[k+1] + VHP[k]


    return forward, surface_value, surface_flux, bottom_flux


@njit
def zonal_velocity_profile(dt2, num_layers, column_depth, dz, dzz, umol, km, uf, bsu, wsu):
    """ 
    Description: Calculates zonal (U) velocity profile
                 Solves for the equation:    dti2 * (KM * U')' - U = -UB
                 
                 dt2 = twice the timestep
                 dz = vertical spacing
                 dzz = staggered vertical spacing
                 km = diffusion of momentum
                 uf = forward step for zonal velocity
                 bsu = zonal bottom stress
                 wsu = zonal wind stress
                 umol = background diffusion coefficient
    
    :return: data array for zonal velocity profile
    """
    A = np.zeros(num_layers)
    C = np.zeros(num_layers)
    VH = np.zeros(num_layers)
    VHP = np.zeros(num_layers)

    A[:-2] = -dt2 * (km[1:-1] + umol) / (dz[:-2] * dzz[:-2] * column_depth * column_depth)
    C[1:-1] = -dt2 * (km[1:-1] + umol) / (dz[1:-1] * dzz[:-2] * column_depth * column_depth)
        

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-dt2 * wsu / (-dz[0] * column_depth) - uf[0]) / (A[0] - 1.)

    for i in range(1, num_layers - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - uf[i]) * VHP[i]

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-dt2 * wsu / (-dz[0] * column_depth) - uf[0]) / (A[0] - 1.)

    CBC = 0.0
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    uf[num_layers - 2] = (C[num_layers - 2] * VHP[num_layers - 3] - uf[num_layers - 2]) / (CBC * dt2 / (-dz[num_layers - 2] * column_depth) - 1. - (VH[num_layers - 3] - 1.) * C[num_layers - 2])

    for i in range(1, num_layers - 1):
        k = num_layers - 1 - i
        uf[k - 1] = VH[k - 1] * uf[k] + VHP[k - 1]

    bsu = -CBC * uf[num_layers - 2]
    
    return uf, bsu
