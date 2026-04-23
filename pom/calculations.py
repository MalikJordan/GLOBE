import numpy as np
from pom.check_diffusion import a, c, vh, vhp

# def density_profile(physical):
    
#     gravity = 9.806
#     vertical_density_profile = np.zeros(physical["water_column"]["num_layers"])

#     pressure = -gravity * 1.025 * physical["vertical_grid"]["dzz"][:-1] * physical["water_column"]["column_depth"] * 0.01

#     cr = 1449.1 + (0.0821*pressure) + (4.55*physical["temperature"]["tb"][:-1]) - (0.045*np.power(physical["temperature"]["tb"][:-1],2)) + (1.34*(physical["salinity"]["sb"][:-1] - 35.0))
#     cr = pressure/np.power(cr,2)
    
#     density = 999.842594 + (6.793952e-02*physical["temperature"]["tb"][:-1]) - (9.095290e-03*np.power(physical["temperature"]["tb"][:-1],2)) \
#             + (1.001685e-04*np.power(physical["temperature"]["tb"][:-1],3)) - (1.120083e-06*np.power(physical["temperature"]["tb"][:-1],4)) + (6.536332e-09*np.power(physical["temperature"]["tb"][:-1],5)) \
#             + (0.824493 - (4.0899e-03*physical["temperature"]["tb"][:-1]) + (7.6438e-05*np.power(physical["temperature"]["tb"][:-1],2))
#                     - (8.2467e-07*np.power(physical["temperature"]["tb"][:-1],3)) + (5.3875e-09*np.power(physical["temperature"]["tb"][:-1],4))) * physical["salinity"]["sb"][:-1] \
#             + (-5.72466e-03 + 1.0227e-04*physical["temperature"]["tb"][:-1] - (1.6546e-06*np.power(physical["temperature"]["tb"][:-1],2))) * (np.power(np.abs(physical["salinity"]["sb"][:-1]),1.5)) \
#                     + (4.8314e-04*np.power(physical["salinity"]["sb"][:-1],2)) + 1.0e05*cr*(1.0 - (2*cr))
    
#     vertical_density_profile[:-1] = (density - 1000.) * 1.e-03
#     vertical_density_profile[-1] = vertical_density_profile[-2]

#     physical["density"] = vertical_density_profile

#     return physical


def density_profile(physical):
    
    gravity = 9.806
    vertical_density_profile = np.zeros(physical["water_column"]["num_layers"])
    pressure = -gravity * 1.025 * physical["vertical_grid"]["dzz"][:-1] * physical["water_column"]["column_depth"] * 0.01

    # density = 999.842594 + (6.793952E-02 * physical["temperature"]["tb"][:-1]) - (9.095290E-03 * np.power(physical["temperature"]["tb"][:-1],2)) \
    #         + (1.001685E-04 * np.power(physical["temperature"]["tb"][:-1],3)) - (1.120083E-06 * np.power(physical["temperature"]["tb"][:-1],4)) \
    #         + (6.536332E-09 * np.power(physical["temperature"]["tb"][:-1],5))
    # density = density + (0.824493 - (4.0899E-03 * physical["temperature"]["tb"][:-1]) + (7.6438E-05 * np.power(physical["temperature"]["tb"][:-1],2)) \
    #         - (8.2467E-07 * np.power(physical["temperature"]["tb"][:-1],3)) + (5.3875E-09 * np.power(physical["temperature"]["tb"][:-1],4))) * physical["salinity"]["sb"][:-1] \
    #         + (-5.72466E-03 + (1.0227E-04 * physical["temperature"]["tb"][:-1]) - (1.6546E-06 * np.power(physical["temperature"]["tb"][:-1],2))) * (np.power(np.abs(physical["salinity"]["sb"][:-1]),1.5)) \
    #         + 4.8314E-4*np.power(physical["salinity"]["sb"][:-1],2)

    density = 999.842594 + (6.793952E-02 * physical["temperature"]["t"][:-1]) - (9.095290E-03 * np.power(physical["temperature"]["t"][:-1],2)) \
            + (1.001685E-04 * np.power(physical["temperature"]["t"][:-1],3)) - (1.120083E-06 * np.power(physical["temperature"]["t"][:-1],4)) \
            + (6.536332E-09 * np.power(physical["temperature"]["t"][:-1],5))
    density = density + (0.824493 - (4.0899E-03 * physical["temperature"]["t"][:-1]) + (7.6438E-05 * np.power(physical["temperature"]["t"][:-1],2)) \
            - (8.2467E-07 * np.power(physical["temperature"]["t"][:-1],3)) + (5.3875E-09 * np.power(physical["temperature"]["t"][:-1],4))) * physical["salinity"]["s"][:-1] \
            + (-5.72466E-03 + (1.0227E-04 * physical["temperature"]["t"][:-1]) - (1.6546E-06 * np.power(physical["temperature"]["t"][:-1],2))) * (np.power(np.abs(physical["salinity"]["s"][:-1]),1.5)) \
            + 4.8314E-4*np.power(physical["salinity"]["s"][:-1],2)

    vertical_density_profile[:-1] = (density - 1000.) * 1.E-03
    vertical_density_profile[-1] = vertical_density_profile[-2]

    physical["density"] = vertical_density_profile

    return physical


def kinetic_energy_profile(physical, pom1d):
    
    # INITIALIZE VARIABLES
    A = np.zeros(physical["water_column"]["num_layers"])
    C = np.zeros(physical["water_column"]["num_layers"])
    VH = np.zeros(physical["water_column"]["num_layers"])
    VHP = np.zeros(physical["water_column"]["num_layers"])

    KN = np.zeros(physical["water_column"]["num_layers"])
    GH = np.zeros(physical["water_column"]["num_layers"])
    SH = np.zeros(physical["water_column"]["num_layers"])
    SM = np.zeros(physical["water_column"]["num_layers"])

    DTEF = np.zeros(physical["water_column"]["num_layers"])
    BPROD = np.zeros(physical["water_column"]["num_layers"])
    PROD = np.zeros(physical["water_column"]["num_layers"])
    SPROD = np.zeros(physical["water_column"]["num_layers"])
    pressure = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   LOCAL ARRAYS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    BOYGR = np.zeros(physical["water_column"]["num_layers"]); CC = np.zeros(physical["water_column"]["num_layers"])
    TEMP1 = np.zeros(physical["water_column"]["num_layers"]); TEMP2 = np.zeros(physical["water_column"]["num_layers"]); TEMP3 = np.zeros(physical["water_column"]["num_layers"])

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
    
    A[1:-1] = -physical["simulation"]["dt2"] * (physical["diffusion"]["kinetic_energy"][2:] + physical["diffusion"]["kinetic_energy"][1:-1] + 2 * pom1d["background_diffusion"]["umol"]) * \
                0.5 / (physical["vertical_grid"]["dzz"][:-2] * physical["vertical_grid"]["dz"][1:-1] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])
    C[1:-1] = -physical["simulation"]["dt2"] * (physical["diffusion"]["kinetic_energy"][:-2] + physical["diffusion"]["kinetic_energy"][1:-1] + 2 * pom1d["background_diffusion"]["umol"]) * \
                0.5 / (physical["vertical_grid"]["dzz"][:-2] * physical["vertical_grid"]["dz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES FOR THE EQUATION
    #   DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    CONST1 = 16.6 ** 0.6666667 * CIWC

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    VH[0] = 0.0
    VHP[0] = np.sqrt(physical["stresses"]["wsu"] ** 2 + physical["stresses"]["wsv"] ** 2) * CONST1
    physical["kinetic_energy"]["kef"][physical["water_column"]["num_layers"] - 1] = 0.5 * np.sqrt((physical["stresses"]["bsu"] + physical["stresses"]["bsu"]) ** 2 +
                                                                (physical["stresses"]["bsv"] + physical["stresses"]["bsv"]) ** 2) * CONST1

    physical["kinetic_energy"]["keb"][1:-1] = np.abs(physical["kinetic_energy"]["keb"][1:-1])
    physical["kinetic_energy"]["kelb"][1:-1] = np.abs(physical["kinetic_energy"]["kelb"][1:-1])
    BOYGR[1:-1] = gravity * (physical["density"][:-2] - physical["density"][1:-1]) / (physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"])
    DTEF[1:-1] = physical["kinetic_energy"]["keb"][1:-1] * np.sqrt(physical["kinetic_energy"]["keb"][1:-1]) / (B1 * physical["kinetic_energy"]["kelb"][1:-1] + SMALL)
    SPROD[1:-1] = .25 * physical["diffusion"]["momentum"][1:-1] * \
                   ((physical["velocity"]["u"][1:-1] + physical["velocity"]["u"][1:-1] - physical["velocity"]["u"][0:-2] - physical["velocity"]["u"][0:-2]) ** 2
                    + (physical["velocity"]["v"][1:-1] + physical["velocity"]["v"][1:-1] - physical["velocity"]["v"][0:-2] - physical["velocity"]["v"][0:-2]) ** 2) / \
                   (physical["vertical_grid"]["dzz"][0:-2] * physical["water_column"]["column_depth"]) ** 2 * CIWC ** 2
    BPROD[1:-1] = physical["diffusion"]["tracers"][1:-1] * BOYGR[1:-1]
    PROD[1:-1] = SPROD[1:-1] + BPROD[1:-1]
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for i in range(1, physical["water_column"]["num_layers"] - 1):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (2. * physical["simulation"]["dt2"] * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (-2. * physical["simulation"]["dt2"] * PROD[i] + C[i] * VHP[i - 1] - physical["kinetic_energy"]["keb"][i]) * VHP[i]
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    for i in range(2, physical["water_column"]["num_layers"]+1):  # 104
        k = physical["water_column"]["num_layers"] - i
        physical["kinetic_energy"]["kef"][k] = VH[k] * physical["kinetic_energy"]["kef"][k + 1] + VHP[k]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SEECTION SOLVES FOR TEH EQUATION
    #   DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   BOUNDARY CONDITIONS
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    VH[0] = 0.
    VHP[0] = 0.
    physical["kinetic_energy"]["kelf"][physical["water_column"]["num_layers"] - 1] = 0.

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP DOWNWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, physical["water_column"]["num_layers"] - 1):
        DTEF[i] = DTEF[i] * (1. + E2 * ((1. / np.abs(physical["vertical_grid"]["z"][i] - physical["vertical_grid"]["z"][0])
                                         + 1. / np.abs(physical["vertical_grid"]["z"][i] - physical["vertical_grid"]["z"][physical["water_column"]["num_layers"]-1]))
                                        * physical["vertical_grid"]["l"][i] / (physical["water_column"]["column_depth"] * von_karman_constant)) ** 2)
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - (physical["simulation"]["dt2"] * DTEF[i] + 1.))
        VH[i] = A[i] * VHP[i]
        VHP[i] = (physical["simulation"]["dt2"] * (- (SPROD[i] + E3 * BPROD[i]) * physical["vertical_grid"]["l"][i] * E1)
                  + C[i] * VHP[i - 1] - physical["kinetic_energy"]["kelb"][i]) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SWEEP UPWARD
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(2, physical["water_column"]["num_layers"]+1):
        k = physical["water_column"]["num_layers"] - i
        physical["kinetic_energy"]["kelf"][k] = VH[k] * physical["kinetic_energy"]["kelf"][k + 1] + VHP[k]

    for i in range(1, physical["water_column"]["num_layers"] - 1):
        if physical["kinetic_energy"]["kef"][i] > SMALL or physical["kinetic_energy"]["kelf"][i] > SMALL:
            continue
        else:
            physical["kinetic_energy"]["kef"][i] = SMALL
            physical["kinetic_energy"]["kelf"][i] = SMALL

    physical["kinetic_energy"]["kef"][:-1] = np.abs(physical["kinetic_energy"]["kef"][:-1])
    physical["kinetic_energy"]["kelf"][:-1] = np.abs(physical["kinetic_energy"]["kelf"][:-1])
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
    physical["vertical_grid"]["l"][0] = 0.
    physical["vertical_grid"]["l"][physical["water_column"]["num_layers"] - 1] = 0.
    GH[0] = 0.
    GH[physical["water_column"]["num_layers"] - 1] = 0.

    physical["vertical_grid"]["l"][1:-1] = physical["kinetic_energy"]["kelf"][1:-1] / physical["kinetic_energy"]["kef"][1:-1]
    GH[1:-1] = np.power(physical["vertical_grid"]["l"][1:-1],2) / physical["kinetic_energy"]["kef"][1:-1] * BOYGR[1:-1]
    for i in range(0, physical["water_column"]["num_layers"]):
        GH[i] = min(GH[i], .028)
    SH = COEF1 / (1. - COEF2 * GH)
    SM = COEF3 + SH * COEF4 * GH
    SM = SM / (1. - COEF5 * GH)

    KN = physical["vertical_grid"]["l"] * np.sqrt(np.abs(physical["kinetic_energy"]["ke"]))
    physical["diffusion"]["kinetic_energy"] = 0.5 * (KN * 0.41 * SM + physical["diffusion"]["kinetic_energy"])
    physical["diffusion"]["momentum"] = 0.5 * (KN * SM + physical["diffusion"]["momentum"])
    physical["diffusion"]["tracers"] = 0.5 * (KN * SH + physical["diffusion"]["tracers"])

    return physical


def meridional_velocity_profile(physical, pom1d):
    """ 
    Description: Calculates meridional (V) velocity profile
                 Solves for the equation:    dti2 * (KM * V')' - V = -VB
    
    :return: data array for meridional velocity profile
    """
    physical["stresses"]["wsv"]
    A = np.zeros(physical["water_column"]["num_layers"])
    C = np.zeros(physical["water_column"]["num_layers"])
    VH = np.zeros(physical["water_column"]["num_layers"])
    VHP = np.zeros(physical["water_column"]["num_layers"])
  
    A[:-2] = -physical["simulation"]["dt2"] * (physical["diffusion"]["momentum"][1:-1] + pom1d["background_diffusion"]["umol"]) / \
                (physical["vertical_grid"]["dz"][:-2] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])
    C[1:-1] = -physical["simulation"]["dt2"] * (physical["diffusion"]["momentum"][1:-1] + pom1d["background_diffusion"]["umol"]) / \
            (physical["vertical_grid"]["dz"][1:-1] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-physical["simulation"]["dt2"] * physical["stresses"]["wsv"] / (-physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - physical["velocity"]["vf"][0]) / (A[0] - 1.)

    # 98 CONTINUE

    for i in range(1, physical["water_column"]["num_layers"] - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - physical["velocity"]["vf"][i]) * VHP[i]

    CBC = 0.0
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   TO RESTORE BOTTOM B.L. DELETE NEXT LINE
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    physical["velocity"]["vf"][physical["water_column"]["num_layers"] - 2] = (C[physical["water_column"]["num_layers"] - 1] * VHP[physical["water_column"]["num_layers"] - 3] - physical["velocity"]["vf"][physical["water_column"]["num_layers"] - 2]) / (
            CBC * physical["simulation"]["dt2"] / (-physical["vertical_grid"]["dz"][physical["water_column"]["num_layers"] - 2] * physical["water_column"]["column_depth"]) - 1. - (VH[physical["water_column"]["num_layers"] - 3] - 1.) * C[physical["water_column"]["num_layers"] - 2])

    for i in range(1, physical["water_column"]["num_layers"] - 1):
        k = physical["water_column"]["num_layers"] - 1 - i
        physical["velocity"]["vf"][k - 1] = VH[k - 1] * physical["velocity"]["vf"][k] + VHP[k - 1]

    physical["stresses"]["bsv"] = -CBC * physical["velocity"]["vf"][physical["water_column"]["num_layers"] - 2]  # 92

    return physical


def mixed_layer_depth(physical):
    small = 1.e-06
    
    for i in range(0,physical["water_column"]["num_layers"]-1):

        if physical["temperature"]["t"][0] > physical["temperature"]["t"][i]+0.2:
            break

        physical["mld"][i] = physical["vertical_grid"]["zz"][i] - (physical["temperature"]["t"][i] + 0.2 - physical["temperature"]["t"][0]) * \
                            (physical["vertical_grid"]["zz"][i] - physical["vertical_grid"]["zz"][i+1]) / \
                               (physical["temperature"]["t"][i] - physical["temperature"]["t"][i+1] + small)
        
    return physical


def temperature_and_salinity_profiles(physical, pom1d, property, case):
    """ 
    Description: Solves for the conservative (temperature & salinity) and non-conservative (BFM state variables) scalars
                 Handles the surface and bottom boundary conditions
    NOTE: Conservative scalars are only calculated when the system is run in prognostic mode
    
    :return: data arrays for the 'property' (temperature, salinity, or BFM state variable)
    """
    # Read property data
    if case == 'Temperature':
        forward = property["tf"]
        surface_value = property["surf"]
        surface_flux = property["surf_flux"]
        bottom_flux = property["bot_flux"]
    elif case == 'Salinity':
        forward = property["sf"]
        surface_value = property["surf"]
        surface_flux = property["surf_flux"]
        bottom_flux = property["bot_flux"]
    elif case == 'BGC':
        forward = property["bf"]
        surface_value = property["surf"]
        surface_flux = property["surf_flux"]
        bottom_flux = property["bot_flux"]
    
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
    if case == 'Temperature':
        nbc  = pom1d["flags"]["nbct"]
        umol = pom1d["background_diffusion"]["umolt"]
    elif case == 'Salinity':
        nbc  = pom1d["flags"]["nbcs"]
        umol = pom1d["background_diffusion"]["umols"]
    elif case == 'BGC':
        nbc  = pom1d["flags"]["nbcbgc"]
        umol = pom1d["background_diffusion"]["umolbgc"]

    A = np.zeros(physical["water_column"]["num_layers"])
    C = np.zeros(physical["water_column"]["num_layers"])
    VH = np.zeros(physical["water_column"]["num_layers"])
    VHP = np.zeros(physical["water_column"]["num_layers"])

    # SW PROFILE
    vertical_radiation_profile = np.zeros(physical["water_column"]["num_layers"])

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
    A[:-2] = -physical["simulation"]["dt2"] * (physical["diffusion"]["tracers"][1:-1] + umol) / (physical["vertical_grid"]["dz"][:-2] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])
    C[1:-1] = -physical["simulation"]["dt2"] * (physical["diffusion"]["tracers"][1:-1] + umol) / (physical["vertical_grid"]["dz"][1:-1] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])

    vertical_radiation_profile[:] = 0.
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   SURFACE BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    # *** PENETRATIVE RADIATION CALCULATION. AT THE BOTTOM ANY UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER.
    
    if nbc == 1:
        VH[0] = A[0] / (A[0] - 1.)
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            VHP[0] = -physical["simulation"]["dt2"] * (surface_flux + 0.) / (-physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - forward[0]
        else:
            VHP[0] = -physical["simulation"]["dt2"] * (surface_flux + physical["swrad"]) / (-physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 2:
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            vertical_radiation_profile[:] = 0. * (RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD1[pom1d["flags"]["ntp"]]) + (1. - RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD2[pom1d["flags"]["ntp"]])))  # ***
        else:
            vertical_radiation_profile[:] = physical["swrad"] * (RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD1[pom1d["flags"]["ntp"]]) + (1. - RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD2[pom1d["flags"]["ntp"]])))  # ***
        vertical_radiation_profile[physical["water_column"]["num_layers"] - 1] = 0.

        VH[0] = A[0] / (A[0] - 1.)
        VHP[0] = physical["simulation"]["dt2"] * (surface_flux + vertical_radiation_profile[0] - vertical_radiation_profile[1]) / (physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - forward[0]
        VHP[0] = VHP[0] / (A[0] - 1.)

    elif nbc == 3:
        VH[0] = 0.
        VHP[0] = surface_value

    elif nbc == 4:
        if case == 'BGC':   # Exclude shortwave radiation in calculations
            vertical_radiation_profile[:] = 0. * (RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD1[pom1d["flags"]["ntp"]]) + (1. - RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD2[pom1d["flags"]["ntp"]])))  # ***
        else:
            vertical_radiation_profile[:] = physical["swrad"] * (RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD1[pom1d["flags"]["ntp"]]) + (1. - RP[pom1d["flags"]["ntp"]] * np.exp(physical["vertical_grid"]["z"][:] * physical["water_column"]["column_depth"] / AD2[pom1d["flags"]["ntp"]])))  # ***
        vertical_radiation_profile[physical["water_column"]["num_layers"] - 1] = 0.

        VH[0] = 0.
        VHP[0] = surface_value

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   THE FOLLOWING SECTION SOLVES THE EQUATION
    #   DT2*(KH*FF')' -FF = -FB
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(1, physical["water_column"]["num_layers"] - 2):
        VHP[i] = 1 / (A[i] + C[i] * (1 - VH[i - 1]) - 1)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - forward[i] + physical["simulation"]["dt2"] * (vertical_radiation_profile[i] - vertical_radiation_profile[i + 1]) / (physical["water_column"]["column_depth"] * physical["vertical_grid"]["dz"][i])) * VHP[i]

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    forward[physical["water_column"]["num_layers"] - 2] = (C[physical["water_column"]["num_layers"] - 2] * VHP[physical["water_column"]["num_layers"] - 3] - forward[physical["water_column"]["num_layers"] - 2] \
                                                                 + (bottom_flux * physical["simulation"]["dt2"] / (physical["vertical_grid"]["dz"][physical["water_column"]["num_layers"] - 2] * physical["water_column"]["column_depth"]))) \
                                                                / (C[physical["water_column"]["num_layers"] - 2] * (1 - VH[physical["water_column"]["num_layers"] - 3]) - 1)
    
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #   APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    for i in range(3, physical["water_column"]["num_layers"]+1):
        k = (physical["water_column"]["num_layers"]) - i
        forward[k] = VH[k] * forward[k+1] + VHP[k]


    # Update property data
    if case == 'Temperature':
        property["tf"] = forward
        property["surf"] = surface_value
        property["surf_flux"] = surface_flux
        property["bot_flux"] = bottom_flux
    elif case == 'Salinity':
        property["sf"] = forward
        property["surf"] = surface_value
        property["surf_flux"] = surface_flux
        property["bot_flux"] = bottom_flux
    elif case == 'BGC':
        property["bf"] = forward
        property["surf"] = surface_value
        property["surf_flux"] = surface_flux
        property["bot_flux"] = bottom_flux


    return property


def zonal_velocity_profile(physical, pom1d):
    """ 
    Description: Calculates zonal (U) velocity profile
                 Solves for the equation:    dti2 * (KM * U')' - U = -UB
    
    :return: data array for zonal velocity profile
    """
    A = np.zeros(physical["water_column"]["num_layers"])
    C = np.zeros(physical["water_column"]["num_layers"])
    VH = np.zeros(physical["water_column"]["num_layers"])
    VHP = np.zeros(physical["water_column"]["num_layers"])

    A[:-2] = -physical["simulation"]["dt2"] * (physical["diffusion"]["momentum"][1:-1] + pom1d["background_diffusion"]["umol"]) / \
                (physical["vertical_grid"]["dz"][:-2] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])
    C[1:-1] = -physical["simulation"]["dt2"] * (physical["diffusion"]["momentum"][1:-1] + pom1d["background_diffusion"]["umol"]) / \
               (physical["vertical_grid"]["dz"][1:-1] * physical["vertical_grid"]["dzz"][:-2] * physical["water_column"]["column_depth"] * physical["water_column"]["column_depth"])
        

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-physical["simulation"]["dt2"] * physical["stresses"]["wsu"] / (-physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - physical["velocity"]["uf"][0]) / (A[0] - 1.)

    for i in range(1, physical["water_column"]["num_layers"] - 2):
        VHP[i] = 1. / (A[i] + C[i] * (1. - VH[i - 1]) - 1.)
        VH[i] = A[i] * VHP[i]
        VHP[i] = (C[i] * VHP[i - 1] - physical["velocity"]["uf"][i]) * VHP[i]

    VH[0] = A[0] / (A[0] - 1.)
    VHP[0] = (-physical["simulation"]["dt2"] * physical["stresses"]["wsu"] / (-physical["vertical_grid"]["dz"][0] * physical["water_column"]["column_depth"]) - physical["velocity"]["uf"][0]) / (A[0] - 1.)

    CBC = 0.0
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    physical["velocity"]["uf"][physical["water_column"]["num_layers"] - 2] = (C[physical["water_column"]["num_layers"] - 2] * VHP[physical["water_column"]["num_layers"] - 3] - physical["velocity"]["uf"][physical["water_column"]["num_layers"] - 2]) / (
            CBC * physical["simulation"]["dt2"] / (-physical["vertical_grid"]["dz"][physical["water_column"]["num_layers"] - 2] * physical["water_column"]["column_depth"]) - 1. - (VH[physical["water_column"]["num_layers"] - 3] - 1.) * C[physical["water_column"]["num_layers"] - 2])
    for i in range(1, physical["water_column"]["num_layers"] - 1):
        k = physical["water_column"]["num_layers"] - 1 - i
        physical["velocity"]["uf"][k - 1] = VH[k - 1] * physical["velocity"]["uf"][k] + VHP[k - 1]
    physical["stresses"]["bsu"] = -CBC * physical["velocity"]["uf"][physical["water_column"]["num_layers"] - 2]  # 92
    
    return physical
