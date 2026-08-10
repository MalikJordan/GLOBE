import os
import sys
import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, nutrient_limitation, temperature_dependence, tracer_elements, switch, calculate_acidity, calculate_Hplus, find_roots_of_f_TA
from fractions import Fraction
np.set_printoptions(precision=20)
class Inorganic():
    """
    
    """

    def __init__(self, abbrev, physical, reactions, **tracer):

        # Variales that will be used later ---------------------------------------------------------------
        num_layers = physical["water_column"]["num_layers"]
        iters = physical["simulation"]["iters"]
        composition = physical["initial_conditions"][abbrev]["composition"]

        # Add important keys ---------------------------------------------------------------
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Add parameters ---------------------------------------------------------------
        # Surface flux
        if "surface_flux" in tracer["parameters"]:
            if tracer["parameters"]["surface_flux"] == True:
                self.surf_flux = 0.

        # CO2 flux
        if "air_sea_flux" in tracer["parameters"]:
            self.air_sea_flux_ids = List.empty_list(types.unicode_type)
            self.air_sea_flux_params = []

            for key,val in tracer["parameters"]["air_sea_flux"].items():
                self.air_sea_flux_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.air_sea_flux_params.append(val)

        # Denitrification
        if "denitrification" in tracer["parameters"]:
            self.denitrification_ids = List.empty_list(types.unicode_type)
            self.denitrification_params = []

            if "convert_o2" in tracer["parameters"]["denitrification"]:
                # If conversion given as a Fraction string (ex: '1/2' instead of 0.5), convert to float
                if isinstance(tracer["parameters"]["denitrification"]["convert_o2"],str):
                    tracer["parameters"]["denitrification"]["convert_o2"] = np.float64(Fraction(tracer["parameters"]["denitrification"]["convert_o2"]))

            for key,val in tracer["parameters"]["denitrification"].items():
                self.denitrification_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.denitrification_params.append(val)

        # Nitrification
        if "nitrification" in tracer["parameters"]:
            self.nitrification_ids = List.empty_list(types.unicode_type)
            self.nitrification_params = []

            if "convert_o2" in tracer["parameters"]["nitrification"]:
                # If conversion given as a Fraction string (ex: '1/2' instead of 0.5), convert to float
                if isinstance(tracer["parameters"]["nitrification"]["convert_o2"],str):
                    tracer["parameters"]["nitrification"]["convert_o2"] = np.float64(Fraction(tracer["parameters"]["nitrification"]["convert_o2"]))

            for key,val in tracer["parameters"]["nitrification"].items():
                self.nitrification_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.nitrification_params.append(val)
        
        # Rearation
        if "reaeration" in tracer["parameters"]:
            self.reaeration_ids = List.empty_list(types.unicode_type)
            self.reaeration_params = []

            for key,val in tracer["parameters"]["reaeration"].items():
                self.reaeration_ids.append(key)
                if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                self.reaeration_params.append(val)

        # Reoxidation
        if "reoxidation" in tracer["parameters"]:
            self.reoxidation_ids = List.empty_list(types.unicode_type)
            self.reoxidation_params = []

            if "convert_o2" in tracer["parameters"]["reoxidation"]:
                # If conversion given as a Fraction string (ex: '1/2' instead of 0.5), convert to float
                if isinstance(tracer["parameters"]["reoxidation"]["convert_o2"],str):
                    tracer["parameters"]["reoxidation"]["convert_o2"] = np.float64(Fraction(tracer["parameters"]["reoxidation"]["convert_o2"]))

            for key,val in tracer["parameters"]["reoxidation"].items():
                self.reoxidation_ids.append(key)
                self.reoxidation_params.append(val)

        # Temperature regulation
        if "temperature_regulation" in tracer["parameters"]:    
            # Initialize temperature regulation factor to array of 1. (if phytoplankton is temperature limited this will be updated later otherwise will stay as 1.)
            if num_layers > 1:  self.temp_regulation_factor = np.ones(num_layers-1, dtype=np.float64)
            else:   self.temp_regulation_factor = np.float64(1.)
            
            if "temp_limited" in tracer["parameters"]["temperature_regulation"]:    self.temp_limited = tracer["parameters"]["temperature_regulation"]["temp_limited"]
            else:   self.temp_limited = False

            if self.temp_limited:
                # Create float option numbers for use in numba typed.List
                if "function" in tracer["parameters"]["temperature_regulation"]:
                    if tracer["parameters"]["temperature_regulation"]["function"] == "arrhenius":   tracer["parameters"]["temperature_regulation"]["function"] = 1
                    elif tracer["parameters"]["temperature_regulation"]["function"] == "eppley":    tracer["parameters"]["temperature_regulation"]["function"] = 2
                    elif tracer["parameters"]["temperature_regulation"]["function"] == "q10":       tracer["parameters"]["temperature_regulation"]["function"] = 3

                # Add temperature regulation factor to dictionary
                # tracer["parameters"]["temperature_regulation"]["temp_regulation_factor"] = np.empty((0,),dtype=np.float64)

                self.temp_reg_ids = List.empty_list(unicode_type)
                self.temp_reg_params = List.empty_list(float64[:])
                
                for key,val in tracer["parameters"]["temperature_regulation"].items():
                    self.temp_reg_ids.append(key)
                    if not isinstance(val, np.ndarray): val = np.array([val],dtype=np.float64)  # Convert type to array of floats for typed.List
                    self.temp_reg_params.append(val)
        else: self.temp_limited = False

        # Add concentrations ---------------------------------------------------------------
        # Concentration array
        self.composition = []
        conc = []
        if len(tracer["composition"]) > 1:    sys.exit("Inorganic: Only one element accepted per inorganic nutrient. Check documentation adn edit input file.")
        elif len(tracer["composition"]) < 1:  sys.exit("Inorganic: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:   pass

        for key in composition: 
            # Add constituent to composition/concentration
            self.composition.append(key)

            # Set initial conditions
            if isinstance(composition[key], str): # Read initial conditions from file
                conc.append( np.fromfile(os.getcwd() + composition[key]) )

            elif isinstance(composition[key], (int,float)): # Create array of initial conditions
                if num_layers == 1: # 0d configuration
                    conc.append(composition[key])
                else: # 1d configuration
                    conc.append(composition[key] * np.ones(num_layers,dtype=np.float64))

            elif isinstance(composition[key], (list,np.ndarray)):   # Already an array of initial conditions
                conc.append(np.array(composition[key]))

        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.conc = np.zeros((len(self.composition),num_layers-1,iters),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const][:-1]
        else:   # Model as single box
            self.conc = np.zeros((len(self.composition),iters),dtype=np.float64)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const]
        self.d_dt = np.zeros_like(self.conc[...,0],dtype=np.float64)
        self.conc_ratio = np.ones_like(self.conc[...,0],dtype=np.float64)
        
        # Add reactions ---------------------------------------------------------------
        self.reactions = []
        for reac in reactions:
            # Add reaction to dictionary
            if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
            else:   consumed = {"empty": "empty"}
            if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
            else:   produced = {"empty": "empty"}
            if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
                self.reactions.append(reac)
   
    
    def inorg(self, bact_limitation_factor, conc, d_dt, tracer_map, z, dz, temperature, salinity, density, wind):
        
        if self.temp_limited:
            self.temp_regulation_factor = temperature_dependence(temperature, self.temp_reg_ids, self.temp_reg_params) # calculate temperature regulation factor for nitrification
        
        for reac in self.reactions:
            if reac["type"] == "co2_flux":      self.co2_flux(self.air_sea_flux_ids, self.air_sea_flux_params, temperature[0], salinity[0], density[0], wind, z[0], conc, d_dt, tracer_map)
            if reac["type"] == "denitrification" and hasattr(self,"denitrification_params"):
                # Calculate hs_limitation factor
                if "hs" in tracer_map:  hs_limitation_factor = nutrient_limitation(conc[tracer_map["hs"][0]], 1.)
                else:   hs_limitation_factor = 1.
                self.denitrification(self.denitrification_ids, self.denitrification_params, self.temp_regulation_factor, hs_limitation_factor, bact_limitation_factor, conc, d_dt, tracer_map)

            if reac["type"] == "nitrification" and hasattr(self,"nitrification_params"):  
                if "o2" in tracer_map: 
                    half_sat = self.nitrification_ids.index("half_sat_oxygen")
                    oxy_limitation_factor = nutrient_limitation(conc[tracer_map["o2"][0]],self.nitrification_params[half_sat])
                else:   oxy_limitation_factor = 1.
                self.nitrification(self.nitrification_ids, self.nitrification_params, self.temp_regulation_factor, oxy_limitation_factor, conc[tracer_map["nh4"][0]], d_dt, tracer_map)
            
            if reac["type"] == "reaeration":    self.surf_flux = self.reaeration(self.reaeration_ids, self.reaeration_params, z, temperature, salinity, wind, conc[tracer_map["o2"][0]], d_dt, tracer_map)
            
            if reac["type"] == "reoxidation" and self.abbrev == "hs":   
                if "o2" in tracer_map: 
                    half_sat = self.reoxidation_ids.index("half_sat_oxygen")
                    oxy_limitation_factor = nutrient_limitation(conc[tracer_map["o2"][0]],self.reoxidation_params[half_sat])
                else:   oxy_limitation_factor = 1.
                self.reoxidation(self.reoxidation_ids, self.reoxidation_params, oxy_limitation_factor, conc[tracer_map["hs"][0]], d_dt, tracer_map)


    @staticmethod
    @njit
    def denitrification(denitrification_ids, denitrification_params, temp_regulation_factor, hs_limitation_factor, bact_limitation_factor, conc, d_dt, tracer_map):

        # Extract parameter indices
        denit_rate = denitrification_ids.index("denitrification_rate")
        # oxic_anoxic_coeff = denitrification_ids.index("oxic_anoxic_coeff")
        # nit_anoxic_coeff = denitrification_ids("nit_anoxic_coeff")
        anoxic_mineralization_rate = denitrification_ids.index("anoxic_mineralization_rate")

        # Get concentrations
        no3 = conc[tracer_map["no3"][0]]
        if "o2" in tracer_map:  o2 = conc[tracer_map["o2"][0]]
        else:   o2 = np.zeros_like(no3)
        if "hs" in tracer_map:  hs = conc[tracer_map["hs"][0]]
        else:   hs = np.zeros_like(no3)
        

        # Convert bacteria limitation
        # bact_limitation_factor /= denitrification_params[oxic_anoxic_coeff]
        bact_limitation_factor /= 0.5   # / 0.5 for oxic-anoxic stoichiometric conversion

        # Calculate denitrification
        denitrification = temp_regulation_factor * hs_limitation_factor * bact_limitation_factor * denitrification_params[denit_rate] * no3 / denitrification_params[anoxic_mineralization_rate]
        
        # Maximum of 0.
        denitrification = np.maximum(np.zeros_like(denitrification),denitrification)

        # Update d_dt
        d_dt[tracer_map["no3"][0]] -= denitrification
        if "n2" in tracer_map:  d_dt[tracer_map["n2"][0]] += denitrification
        if "hs" in tracer_map:
            # x = -(o2 - hs) / denitrification_params[oxic_anoxic_coeff]
            x = -(o2 - hs) / 0.5    # / 0.5 for oxic-anoxic stoichiometric conversion
            y = switch(x)
            # convert = denitrification_params[oxic_anoxic_coeff] * denitrification_params[nit_anoxic_coeff] * y
            convert = 0.5 * 1.25 * y    # * 0.5 for oxic-anoxic stoichiometric conversion, * 1.25 for nit-anoxic stoichiometric conversion
            d_dt[tracer_map["hs"][0]] -= convert * denitrification

    
    @staticmethod
    @njit
    def nitrification(nitrification_ids, nitrification_params, temp_regulation_factor, oxy_limitation_factor, nh4, d_dt, tracer_map):
        """
        Definition:: Calculates nitrification rate
        """
        # Extract parameter indices
        nit_rate = nitrification_ids.index("nitrification_rate")
        convert = False
        if "convert_o2" in nitrification_ids:
            convert = True
            convert_o2 = nitrification_ids.index("convert_o2")

        # Calculate nitrification
        nitrification = temp_regulation_factor * oxy_limitation_factor * nitrification_params[nit_rate] * nh4
        
        # Maximum of 0.
        nitrification = np.maximum(np.zeros_like(nitrification),nitrification)

        d_dt[tracer_map["nh4"][0]] -= nitrification
        d_dt[tracer_map["no3"][0]] += nitrification

        if "o2" in tracer_map:
            if convert: d_dt[tracer_map["o2"][0]] -= nitrification * nitrification_params[convert_o2]
            else:       d_dt[tracer_map["o2"][0]] -= nitrification


    @staticmethod
    @njit
    def reaeration(reaeration_ids, reaeration_params, dz, temperature, salinity, wind, o2, d_dt, tracer_map):

        # Extract parameter indices
        d = reaeration_ids.index("d")
        k1 = reaeration_ids.index("k1")
        k2 = reaeration_ids.index("k2")
        k3 = reaeration_ids.index("k3")
        k4 = reaeration_ids.index("k4")
        schmidt = reaeration_ids.index("schmidt")
    
        # Calculate absolute temperature divided by 100
        abt = (temperature + 273.15) / 100.

        # Calculate theoretical oxygen saturatino for temp + salt and conver into proper units of [mmol O2 / m3]
        oxy_sat = np.exp(-173.4292 + (249.6339/abt) + (143.3483*np.log(abt))-(21.8492*abt) + salinity*(-0.033096 + 0.014259*abt - 0.0017*(abt**2)))/(24.4665E-3)

        # Use this one with BFM17 0D
        # oxy_sat = np.exp(-173.4292 + (249.6339/abt) + (143.3483*np.log(abt))-(21.8492*abt) + salinity*(-0.033096 + 0.014259*abt - 0.0017*(abt**2)))*44.661

        # Calculate Schmidt number, ratio between the kinematic viscosity and the molecular diffusivity of CO2
        schmidt_number = reaeration_params[k1] - ( reaeration_params[k2]*temperature ) + ( reaeration_params[k3]*(temperature**2) ) - ( reaeration_params[k4]*(temperature**3) )
        schmidt_ratio = reaeration_params[schmidt] / schmidt_number

        # Schmidt ratio limited to 0
        schmidt_ratio = np.maximum(np.zeros_like(schmidt_ratio),schmidt_ratio)

        # Calculate wind dependency, including conversion from cm/hr to m/s
        wind_dependency = ( reaeration_params[d]*(wind**2))*np.sqrt(schmidt_ratio)

        # Convert from cm/hr to m/day
        wind_dependency = wind_dependency * 0.01 * 24

        # Calculate flux of o2
        d_o2 = wind_dependency * (oxy_sat - o2)     # total

        surf_o2 = np.zeros_like(dz)
        surf_o2[0] = d_o2[0] / dz[0]    # average over depth
        
        # Update d_dt
        d_dt[tracer_map["o2"][0]] += surf_o2

        return d_o2[0]
    

    @staticmethod
    @njit
    def reoxidation(reoxidation_ids, reoxidation_params, oxy_limitation_factor, hs, d_dt, tracer_map):
        """
        Definition:: Calculation reoxidation rate of reduction equivalents
        """
        # Extract parameter indices
        reox_rate = reoxidation_ids.index("reoxidation_rate")
        convert = False
        if "convert_o2" in reoxidation_ids:
            convert = True
            convert_o2 = reoxidation_ids.index("convert_o2")

        reoxidation = reoxidation_params[reox_rate] * oxy_limitation_factor * hs

        # Update d_dt
        d_dt[tracer_map["hs"][0]] -= reoxidation
        if "o2" in tracer_map:
            if convert: d_dt[tracer_map["o2"][0]] -= reoxidation * reoxidation_params[convert_o2]
            else:       d_dt[tracer_map["o2"][0]] -= reoxidation


    def caco3_saturation():
        pass
    
    
    @staticmethod
    # @njit
    def co2_flux(air_sea_flux_ids, air_sea_flux_params, temperature, salinity, density, wind, dz, conc, d_dt, tracer_map):
        """
        Definition:: Calculates pH and rate of co2 air-sea flux
        """

        # Extract parameter indices
        d = air_sea_flux_ids.index("d")
        c1 = air_sea_flux_ids.index("c1")
        c2 = air_sea_flux_ids.index("c2")
        c3 = air_sea_flux_ids.index("c3")
        c4 = air_sea_flux_ids.index("c4")
        pco2_air = air_sea_flux_ids.index("atmospheric_pco2")
        schmidt = air_sea_flux_ids.index("schmidt")
        pH_previous = air_sea_flux_ids.index("pH")
        max_iters = air_sea_flux_ids.index("OCMIP_max_iters")
        max_delta_pH = air_sea_flux_ids.index("OCMIP_max_delta_pH")
        accuracy = air_sea_flux_ids.index("OCMIP_accuracy")

        # Get concentrations
        dic = conc[tracer_map["co2"][0]][0]
        ta = conc[tracer_map["ta"][0]][0]

        # Calculate Schmidt number, ratio between the kinematic viscosity and the molecular diffusivity of carbon dioxide
        # schmidt_number = np.float64((air_sea_flux_params[c1] - air_sea_flux_params[c2]*temperature + air_sea_flux_params[c3]*(temperature**2) - air_sea_flux_params[c4]*(temperature**3)))
        schmidt_number = (air_sea_flux_params[c1] - air_sea_flux_params[c2]*temperature + air_sea_flux_params[c3]*(temperature**2) - air_sea_flux_params[c4]*(temperature**3))
        
        # Schmidt_ratio is limited to 0 when T > 40 °C
        # sr = np.float64(air_sea_flux_params[schmidt]/schmidt_number)
        sr = air_sea_flux_params[schmidt][0]/schmidt_number[0]
        schmidt_ratio = max(0., sr)

        # Calculate solubility and acidity constants
        k0, k1, k1p, k2, k2p, k3p, ksi, kw, ks, kf, kb, ken, bt, st, ft, pt, sit, ldic, alk = calculate_acidity(temperature, salinity, density, wind, air_sea_flux_params[d], schmidt_ratio, dic, ta)

        # calculate [H+] total when DIC and TA are known
        small_interval = (air_sea_flux_params[pH_previous]>4.0 and air_sea_flux_params[pH_previous]<9.0)
        if small_interval:
            h1 = 10.0**(-(air_sea_flux_params[pH_previous] + air_sea_flux_params[max_delta_pH]))
            h2 = 10.0**(-(air_sea_flux_params[pH_previous] - air_sea_flux_params[max_delta_pH]))
            h_plus, error = find_roots_of_f_TA(h1, h2, air_sea_flux_params[accuracy], air_sea_flux_params[max_iters], k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
        if not small_interval or error>0:
            h1 = 10.0**(-11.0)
            h2 = 10.0**(-2.0)
            h_plus, error = find_roots_of_f_TA(h1, h2, air_sea_flux_params[accuracy], air_sea_flux_params[max_iters], k1, k2, k1p, k2p, k3p, ksi, kw, ks, kf, kb, bt, st, ft, pt, sit, ldic, alk)
    
        # Derive [co2] and Compute other diagnostic variables (hco3, co3 and ph) and pco2, the co2 partial pressure in the water
        h_plus2 = h_plus*h_plus
        co2 = ldic*h_plus2/(h_plus2 + k1*h_plus + k1*k2)
        pco2_sea = co2/k0
        pH = -np.log10(h_plus)
        hco3 = k1*co2/h_plus
        co3 = k2*hco3/h_plus
            
        # Convert partial pressure of oceanic CO2  from atm to uatm (1.E06)
        pco2_sea = pco2_sea * 1.E06
        
        # Flux co2 in mmol C m^-2 s^-1
        air_sea_flux = ken*(air_sea_flux_params[pco2_air] - pco2_sea)*k0*density/1000.0
        
        # Convert flux to units of mg C m^-3 s^-1
        # 12 = stoichiometric coefficient between Carbon and O2 [C:O2 = 1:12]
        air_sea_flux = air_sea_flux*12/dz

        # # Zero non-surface layers
        # if len(do3cdt_air_sea_flux) > 1:    do3cdt_air_sea_flux[1:] = 0.

        # Update pH in list
        air_sea_flux_params[pH_previous] = pH

        # Create co2 flux array
        co2_flux = np.zeros_like(conc[tracer_map["co2"][0]])
        co2_flux[0] = air_sea_flux

        # Update d_dt
        d_dt[tracer_map["co2"][0]] += co2_flux

        return air_sea_flux[0]