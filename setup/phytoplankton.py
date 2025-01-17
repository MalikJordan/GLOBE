import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutrient_limitation, temperature_regulation, tracer_elements

class Phytoplankton():
    """
    
    """

    def __init__(self, abbrev, iters, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        
        # Nutrient limitation
        self.nutrients = []
        self.nutrient_half_sat = []
        self.nutrient_limitation_type = tracer["parameters"]["nutrient_limitation"]
        self.nutrient_limitation_factor = 0.
        self.nutrient_limitation = tracer["parameters"]["nutrient_limitation"]
        # for key in tracer["parameters"]["nutrient_limitation"]:
        #     self.nutrient_limitation[key] = tracer["parameters"]["nutrient_limitation"][key]

        # Temperature regulation
        self.temp_limited = tracer["parameters"]["temp_limited"]
        if self.temp_limited:
            self.q10 = tracer["parameters"]["q10"]
        self.fT = 1.

        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Phytoplankton: Element required for " + self.name + ". Check documentation adn edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    self.composition.append(key)
                    conc.append(tracer["composition"][key])
                else:
                    sys.exit("Phytoplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
        
        hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        for i in range(0,len(conc)):
            hold[i,0] = np.array(conc[i])
        self.conc = np.array(hold)
        self.d_dt = np.zeros_like(conc)
        self.conc_ratio = np.zeros_like(conc)

        # Production 
        # self.exu = np.ones_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.exu = np.zeros_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.gpp = np.ones_like(self.conc[0,...],dtype=float)   # Gross Primary Production (Initialized to 1 for use in exudation and respiration)
        self.lys = np.zeros_like(self.conc[0,...],dtype=float)  # Lysis (carbon)
        self.npp = np.zeros_like(self.conc[0,...],dtype=float)  # Net Primary Production
        self.rsp = np.zeros_like(self.conc[0,...],dtype=float)  # Respiration

        # Add relevant reactions
        self.reactions = []
        for reac in reactions:
            # Add reaction to dictionary
            if "consumed" in reac and reac["consumed"] != None:    consumed = reac["consumed"]
            else:   consumed = {"empty": "empty"}
            if "produced" in reac and reac["produced"] != None:    produced = reac["produced"]
            else:   produced = {"empty": "empty"}
            if ( abbrev in consumed.keys() ) or ( abbrev in produced.keys() ):
                self.reactions.append(reac)

        # Reorder reactions
        # uptake / gpp --> exu --> mortality --> respiration --> synthesis
        self.reactions = [item for item in self.reactions if item["type"] == "respiration"] + [item for item in self.reactions if item["type"] != "respiration"]
        self.reactions = [item for item in self.reactions if item["type"] == "mortality"] + [item for item in self.reactions if item["type"] != "mortality"]
        self.reactions = [item for item in self.reactions if item["type"] == "exudation"] + [item for item in self.reactions if item["type"] != "exudation"]
        self.reactions = [item for item in self.reactions if item["type"] == "gross_primary_production"] + [item for item in self.reactions if item["type"] != "gross_primary_production"]
        self.reactions = [item for item in self.reactions if item["type"] == "uptake"] + [item for item in self.reactions if item["type"] != "uptake"]

        # Switch to determine if it is necessary to calculate growth parameters
        self.growth_switch = False
        for reac in self.reactions:
            if reac["type"] == "gross_primary_production" or reac["type"] == "uptake":
                self.growth_switch = True
                break

    def phyto(self, iter, base_element, base_temp, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, tracers):
        
        # Reset production variables
        # self.exu = np.ones_like((self.conc[0,...,0]))   # Exudation (Initialzied to 1 for use in respiration)
        # self.gpp = np.ones_like(self.conc[0,...],dtype=float)   # Gross Primary Production
        # self.lys = np.zeros_like(self.conc[0,...,0])  # Lysis (carbon)
        # self.npp = np.zeros_like(self.conc[0,...,0])  # Net Primary Production
        # self.rsp = np.zeros_like(self.conc[0,...,0])  # Respiration

        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(base_element, iter, tracers)

        # Calculate temp regulation factor        
        if self.temp_limited:
            self.fT = temperature_regulation(base_temp, temperature, self.q10)

        # Calculate growth parameters
        # if self.growth_switch:

        
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "chlorophyll_synthesis":         self.chlorophyll_synthesis(iter, reac["parameters"], coordinates, surface_PAR)
            if reac["type"] == "exudation":                     self.exudation(iter, reac["parameters"], c, p, tracers)
            if reac["type"] == "gross_primary_production":      irrad = self.gross_primary_production(iter, reac["parameters"], coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, tracers)
            if reac["type"] == "mortality":                     self.mortality(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "respiration":                   self.respiration(iter, reac["parameters"], c, p, tracers)
            if reac["type"] == "uptake":                        fI = self.uptake(iter, reac["parameters"], coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, ec, ep, ic, ip, tracers)

        # Calculate net primary production
        self.npp[iter] = np.maximum(np.zeros_like(self.gpp[iter]), (self.gpp[iter] - self.exu[iter] - self.lys[iter] - self.rsp[iter]))
        if iter == 240:
            avg_gpp = np.average(self.gpp[0:iter])
            avg_exu = np.average(self.exu[0:iter])
            x=1

        # return fI, irrad

    def add_nutrient(self, nutrient, half_sat):
        "Nitrate-Phosphate (N-P) co-limitation"
        if nutrient == "no3" or nutrient == "po4":
            self.nutrients.append(nutrient)
            self.nutrient_half_sat.append(half_sat)


    def calculate_nutrient_limitation(self, base_element, iter, tracers):
        """
        Definition:: Calculates nutrient limitation factor as either a minimum, product, or sum of all nutrients which limit phytoplankton growth
        """
        # nutrient_conc = np.zeros(len(self.nutrients),dtype=np.ndarray)
        # for key in tracers:
        #     if key in self.nutrients:
        #         i = self.nutrients.index(key)
        #         nutrient_conc[i] = np.array(tracers[key].conc[...,iter])
        # lim = nutrient_limitation(nutrient_conc, self.nutrient_half_sat)

        # if len(lim) > 0:
        #     if self.nutrient_limitation_type == "minimum":
        #         self.nutrient_limitation_factor = np.min(lim, axis=0)
        #     elif self.nutrient_limitation_type == "product":
        #         self.nutrient_limitation_factor = np.prod(lim, axis=0)
        #     elif self.nutrient_limitation_type == "sum":
        #         self.nutrient_limitation_factor == np.sum(lim,  axis=0)


        fN = np.zeros(len(self.nutrient_limitation["nutrients"]),dtype=np.ndarray)

        if self.nutrient_limitation["type"] == "internal":
            # Get index of base element
            index = self.composition.index(base_element)

            # Calculate concentration ratios
            concentration_ratio(iter, index, self)
            for nut in self.nutrient_limitation["nutrients"]:
                i = self.nutrient_limitation["nutrients"].index(nut)
                # element = str(tracers[nut].composition[0])
                element = tracers[nut].composition[0]
                element_index = self.composition.index(element)

                func = ( self.conc_ratio[element_index] - self.nutrient_limitation["min_quota"][i]) / ( self.nutrient_limitation["opt_quota"][i] - self.nutrient_limitation["min_quota"][i] )
                func = np.maximum(1.E-20*np.ones_like(func), func)
                fN[i] = np.minimum(np.ones_like(func), func)

        elif self.nutrient_limitation["type"] == "external":
            for nut in self.nutrient_limitation["nutrients"]:
                i = self.nutrient_limitation["nutrients"].index(nut)
                # element = str(tracers[nut].composition[0])
                element = tracers[nut].composition[0]
                element_index = self.composition.index(element)

                func = tracers[nut].conc[...,iter] / ( self.nutrient_half_sat[i] + tracers[nut].conc[...,iter] )
                fN[i] = np.maximum(1.E-20*np.ones_like(func), func)
        else:
            sys.exit("Nutrient limitation type not recognized. Check documentation and edit input file.")
        
        if "colimitation" in self.nutrient_limitation.keys():   # Multiple nutrients available
            if self.nutrient_limitation["colimitation"] == "minimum":
                self.nutrient_limitation_factor = np.min(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "product":
                self.nutrient_limitation_factor = np.prod(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "sum":
                self.nutrient_limitation_factor == np.sum(fN,  axis=0)
            else:
                sys.exit("Nutrient colimitation not recognized. Check documentation and edit input file.")
        else:   # Only one nutrient available
            self.nutrient_limitation_factor = fN


    def chlorophyll_synthesis(self, iter, parameters, coordinates, surface_PAR):
        """
        Definition:: Calculates chlorophyll synthesis based on gains and losses of phytoplankton carbon
        """
        # Locate indexes of carbon and chlorophyll constituents
        carbon_index = self.composition.index("c")
        chl_index = self.composition.index("chl")

        # Get carbon and chlorophyll concentrations
        phyto_c = self.conc[carbon_index][iter]
        phyto_chl = self.conc[chl_index][iter]

        # Calculate chlorophyll-carbon ratio
        pl_pc = phyto_chl / phyto_c

        # Calculate irradiance
        k_PAR = light_attenuation(parameters, phyto_c)
        irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)

        # Calculate chlorophyll regulation (rho_chl)
        factor = ( 1 - parameters["activity_respiration_frac"] ) / ( parameters["initial_PI_slope"] * irrad * phyto_chl + 1E-20 )    # 1E-20 to prevent divide by zero error
        rho_chl = parameters["max_chl_c"] * np.minimum(np.ones_like(phyto_chl), np.maximum(np.zeros_like(phyto_chl),factor * (self.gpp[iter] - self.exu[iter])))
        # rho_chl = parameters["max_chl_c"] * np.minimum(np.ones_like(phyto_chl), factor * (self.gpp[iter] - self.exu[iter]))

        # Calculate chlorophyll synthesis
        synthesis = rho_chl * ( 1 - parameters["activity_respiration_frac"] ) * ( self.gpp[iter] - self.exu[iter] ) - pl_pc * ( self.lys[iter] + self.rsp[iter] )
        # synthesis = np.maximum(np.zeros_like(phyto_chl),synthesis)  # synthesis >= 0

        # Update d_dt
        self.d_dt[chl_index] += synthesis


    
    def exudation(self, iter, parameters, c, p, tracers):
        """
        Definition:: Calculates the amount of unassimilated carbon released by phytoplankton into dissolved carbon pool
        """
        # Locate index of carbon constituent 
        carbon_index_phyto = self.composition.index("c")
        carbon_index_om = tracers[p[0]].composition.index("c")

        # Exudation
        # exudation = ( parameters["excreted_fraction"] + ( 1 - parameters["excreted_fraction"] ) * ( 1 - self.nutrient_limitation_factor ) ) * self.gpp[iter]
        exudation = ( 1 - parameters["excreted_fraction"] ) * ( 1 - self.nutrient_limitation_factor ) * self.gpp[iter]

        # Update d_dt
        tracers[c[0]].d_dt[carbon_index_phyto] -= exudation
        tracers[p[0]].d_dt[carbon_index_om] += exudation

        # Update exu variable 
        self.exu[iter] = exudation
    

    def gross_primary_production(self, iter, parameters, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, tracers):
        """
        Definition:: Calculate gross primary production
        """
        
        # Get carbon concentration
        carbon_index = self.composition.index("c")
        phytoc = self.conc[carbon_index][iter]
        
        # Locate index of chlorophyll constituent if present
        if "chl" in self.composition:   
            chl_index = self.composition.index("chl")
            phytol = self.conc[chl_index][iter]
            pl_pc = phytol / phytoc     # Chl:C ratio (used in light limitation)
        else:
            pl_pc = 1.

        # Calculate growth rate
        if parameters["max_photo_rate"] == "variable":
            Vm = max_growth_rate(parameters, temperature)
        else:
            Vm = parameters["max_photo_rate"]

        # Calculate light limitation
        # if parameters["light_limitation"] == "variable":
        if parameters["light_limitation"] in ["monod", "platt", "smith"]:
            # k_PAR = light_attenuation(parameters, phytoc)
            k_PAR = light_attenuation(parameters, phytol)
            irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)
            # fI = light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm)
            fI = light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm)
            # fI = light_limitation(parameters, coordinates, irrad, k_PAR, pl_pc, Vm)
        else:
            fI = parameters["light_limitation"]

        if fI == 1:
            x=1
        # Calculate gross primary production
        gpp = Vm * self.fT * fI * phytoc

        # Update d_dt
        self.d_dt[carbon_index] += gpp
        if "o2" in p:   tracers["o2"].d_dt += gpp / parameters["mw_carbon"]
        if "co2"in c:   tracers["co2"].d_dt -= gpp

        # Update gpp variable
        self.gpp[iter] = gpp

        return irrad


    def lysis(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculate lysis rate of phytoplankton to organic matter pools
        """
        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        tc = np.array(tracers[c].conc[ic][iter])
        tp = np.array(tracers[p].conc[ip][iter])
        
        # Calculate carbon ratios
        carbon_index = self.composition.index('c')
        if 'n' in self.composition:
            n = self.composition.index('n')
            nc = self.conc[n] / self.conc[carbon_index]
        else:
            nc = 1.
        if 'p' in self.composition: 
            p = self.composition.index('p')
            pc = self.conc[p] / self.conc[carbon_index]
        else:
            pc = 1.
        
        # Calculate fraction of lysis released to dissolved pool
        min_quota = np.minimum(parameters["min_nitrogen_quota"]/nc, parameters["min_phosphorus_quota"]/pc)
        factor = np.minimum(np.ones_like(min_quota), min_quota)

        # Calculate nutrient stress limitation
        lim = nutrient_limitation(parameters["nutrient_stress_threshold"], self.nutrient_limitation_factor)

        # Calculate lysis rate
        if parameters["type"] == "particulate": lysis = ( 1 - factor ) * ( lim * parameters["max_lysis_rate"] )
        elif parameters["type"] == "dissolved": lysis = factor * ( lim * parameters["max_lysis_rate"] )

        # Calculate concentration ratios
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        tracers[c].d_dt -= ec * tracers[c].conc_ratio * lysis
        tracers[p].d_dt += ep * tracers[p].conc_ratio * lysis

        # Update lys variable
        self.lys += tracers[c].conc_ratio[carbon_index] * lysis


    def mortality(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates the non-grazing mortality of planktoninc species
                     Lysis is parameterized as a quadratic mortality rate (the carbon constituent of lysis is added to the 'lys' variable for the calculation of net primary production)
        """
        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        tc = np.array(tracers[c].conc[ic][iter])

        # Calculate mortality rate
        mortality = ( parameters["mortality_rate"][0] * tc ) + ( parameters["mortality_rate"][1] * (tc**2) )

        # Temperature regulation
        if self.temp_limited:
            mortality = mortality * self.fT

        # Oxygen limitation
        if "oxygen_limited" in parameters and parameters["oxygen_limited"]:
            fO = nutrient_limitation(tc,parameters["half_sat_oxygen"])
            mortality = mortality * (1 - fO)

        # Calculate concentration ratios
        concentration_ratio(iter, ic, tracers[c])
        concentration_ratio(iter, ip, tracers[p])

        # Update d_dt
        if "partition" in parameters:   # Mortality can be partitioned between dissolved and particulate detrital pools
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
            tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality * parameters["partition"]

        else:
            tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
            tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality

        # Update 'lys' variable if phytoplankton contains a carbon constituent
        if "c" in self.composition:
            carbon_index = self.composition.index("c")
            self.lys[iter] += parameters["mortality_rate"][1] * (self.conc[carbon_index][iter]**2)
        

    def respiration(self, iter, parameters, c, p, tracers):
        """
        Definition:: Calculates phytoplankton respiration
        """
        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        phyto = self.conc[carbon_index][iter]

        # Respiration
        respiration = parameters["respiration_rate"] * self.fT * phyto + parameters["activity_respiration_frac"] * (self.gpp[iter] - self.exu[iter])

        # Update d_dt
        self.d_dt[carbon_index] -= respiration
        if "o2" in c:   tracers["o2"].d_dt -= respiration / parameters["mw_carbon"]
        if "co2" in p:  tracers["co2"].d_dt += respiration

        # Update rsp variable
        self.rsp[iter] = respiration


    def uptake(self, iter, parameters, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, ec, ep, ic, ip, tracers):
        
        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]
        tc = np.array(tracers[c].conc[ic][iter])

        # Get concentration of nutrient in phytoplankton
        nutrient_index = list(ep).index(1.)
        tp = np.array(tracers[p].conc[nutrient_index][iter])

        # Get carbon concentration
        if "c" in self.composition:
            carbon_index = self.composition.index("c")
            phytoc = self.conc[carbon_index][iter]
        else:
            phytoc = np.ones_like(tp)
        
        # Locate index of chlorophyll constituent if present
        if "chl" in self.composition:   
            chl_index = self.composition.index("chl")
            phyto = self.conc[chl_index][iter]
            pl_pc = phyto / phytoc     # Chl:C ratio (used in light limitation)
        else:
            phyto = tp
            pl_pc = 1.

        # Calculate growth rate
        if parameters["max_photo_rate"] == "variable":
            Vm = max_growth_rate(parameters, temperature)
        else:
            Vm = parameters["max_photo_rate"]
        
        # Calculate light limitation
        # if parameters["light_limitation"] == "variable":
        if parameters["light_limitation"] in ["monod", "platt", "smith"]:
            # k_PAR = light_attenuation(parameters, phytoc)
            k_PAR = light_attenuation(parameters, phyto)
            irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)
            # fI = light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm)
            fI = light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm)
        else:
            fI = parameters["light_limitation"]

        # Calculate nutrient limitation
        if parameters["nutrient_limitation"] == "variable":
            fN = nutrient_limitation(tc, parameters["half_sat_nutrient"])
        else:
            fN = parameters["nutrient_limitation"]
        
        # Calculate nutrient uptake
        uptake = Vm * self.fT * fN * fI * tp

        # Calculate concentration ratio
        # concentration_ratio(iter, ic, tracers[c])
        # concentration_ratio(iter, ip, tracers[p])
        
        # Update d_dt
        tracers[c].d_dt -= np.array(ec) * uptake
        tracers[p].d_dt += ep * uptake
        # tracers[c].d_dt -= ec * tracers[c].conc_ratio * uptake
        # tracers[p].d_dt += ep * tracers[p].conc_ratio * uptake
        
        return fI