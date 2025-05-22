import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutrient_limitation, temperature_regulation, tracer_elements, switch

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

        # Intracellular nutrinet quotas
        if "cell_quota" in tracer["parameters"]:
            self.cell_quota = tracer["parameters"]["cell_quota"]

        # Light limitation
        if "light_attenuation" in tracer["parameters"]:
            self.light_attenuation = tracer["parameters"]["light_attenuation"]
        else:
            self.light_attenuation = 0.

        # Temperature regulation
        self.temp_limited = tracer["parameters"]["temp_limited"]
        if self.temp_limited:
            self.q10 = tracer["parameters"]["q10"]
        self.temp_regulation_factor = 1.

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
        # self.conc_ratio = np.zeros_like(conc)
        self.conc_ratio = np.copy(self.cell_quota["opt"])

        # Production 
        # self.exu = np.ones_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.exu = np.zeros_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.gpp = np.zeros_like(self.conc[0,...],dtype=float)   # Gross Primary Production (Initialized to 1 for use in exudation and respiration)
        self.lys = np.zeros_like(self.conc[0,...],dtype=float)  # Lysis (carbon)
        self.npp = np.zeros_like(self.conc[0,...],dtype=float)  # Net Primary Production
        self.psn = np.zeros_like(self.conc[0,...],dtype=float)  # Photosynthesis
        self.rsp = np.zeros_like(self.conc[0,...],dtype=float)  # Respiration

        self.uptn = np.zeros_like(self.conc[0,...],dtype=float) # Nitrogen uptake
        self.uptp = np.zeros_like(self.conc[0,...],dtype=float) # Phosophorus uptake

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
        self.reactions = [item for item in self.reactions if item["type"] == "photosynthesis"] + [item for item in self.reactions if item["type"] != "photosynthesis"]

        # Switch to determine if it is necessary to calculate growth parameters
        self.growth_switch = False
        for reac in self.reactions:
            if reac["type"] == "gross_primary_production" or reac["type"] == "uptake":
                self.growth_switch = True
                break

        # Boolean to determine whether activity and basal respiration will be calculated for chlorophyll synthesis
        self.calc_respiration = False
        for reac in self.reactions:
            if reac["type"] == "respiration":
                self.calc_respiration = True
                break

    def phyto(self, iter, base_element, base_temp, light_attenuation_water, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, tracers):
        
        # Reset production variables
        # self.exu = np.zeros_like((self.conc[0,...,0]))   # Exudation (Initialzied to 1 for use in respiration)
        # self.gpp = np.zeros_like(self.conc[0,...],dtype=float)   # Gross Primary Production
        # self.lys = np.zeros_like(self.conc[0,...,0])  # Lysis (carbon)
        # self.npp = np.zeros_like(self.conc[0,...,0])  # Net Primary Production
        # self.psn = np.zeros_like(self.conc[0,...,0])  # Photosynthesis
        # self.rsp = np.zeros_like(self.conc[0,...,0])  # Respiration

        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(base_element, iter, tracers)

        # Calculate temp regulation factor        
        if self.temp_limited:
            self.temp_regulation_factor = temperature_regulation(base_temp, temperature, self.q10)

        check_conc = self.conc[:,iter]

        # Calculate light extiction coefficient
        k_PAR = light_attenuation(self.abbrev, iter, base_element, light_attenuation_water, tracers)
        
        if iter % 50 == 0:
            x=1
        x=1
        # Calculate rates required for net primary production
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "exudation":                     self.exudation(iter, reac["parameters"], c, p, tracers)
            if reac["type"] == "gross_primary_production":      self.gross_primary_production(iter, reac["parameters"], c, p, tracers)
            if reac["type"] == "lysis":                         self.lysis(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "photosynthesis":                fI, irr = self.photosynthesis(iter, reac["parameters"], base_element, coordinates, dz, k_PAR, temperature, surface_PAR, tracers)
            if reac["type"] == "respiration":                   activity_respiration, basal_respiration = self.respiration(iter, reac["parameters"], c, p, tracers)

        # Calculate net primary production
        self.net_primary_production(iter, tracers)

        # Calculate remaining rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "chlorophyll_synthesis":         
                if self.calc_respiration:   self.chlorophyll_synthesis(iter, reac["parameters"], activity_respiration, basal_respiration, coordinates, irr, k_PAR, surface_PAR, tracers)
                else:                       self.chlorophyll_synthesis(iter, reac["parameters"], 0., 0., coordinates, irr, k_PAR, surface_PAR, tracers)
            if reac["type"] == "uptake":                        self.uptake(iter, reac["parameters"], c, p, ec, ep, ic, tracers)
            
        
        # for reac in self.reactions:
        #     c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
        #     if reac["type"] == "chlorophyll_synthesis":         self.chlorophyll_synthesis(iter, reac["parameters"], coordinates, surface_PAR)
        #     if reac["type"] == "exudation":                     self.exudation(iter, reac["parameters"], c, p, tracers)
        #     if reac["type"] == "gross_primary_production":      irrad = self.gross_primary_production(iter, reac["parameters"], coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, tracers)
        #     if reac["type"] == "mortality":                     self.mortality(iter, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
        #     if reac["type"] == "respiration":                   self.respiration(iter, reac["parameters"], c, p, tracers)
        #     if reac["type"] == "uptake":                        fI = self.uptake(iter, reac["parameters"], coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, ec, ep, ic, ip, tracers)

        # Calculate net primary production
        # self.npp[iter] = np.maximum(np.zeros_like(self.gpp[iter]), (self.gpp[iter] - self.exu[iter] - self.lys[iter] - self.rsp[iter]))
        # if iter == 240:
        #     avg_gpp = np.average(self.gpp[0:iter])
        #     avg_exu = np.average(self.exu[0:iter])
        #     x=1

        # return fI, irrad, self.nutrient_limitation_factor

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
            # concentration_ratio(iter, index, self)

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


    def chlorophyll_synthesis(self, iter, parameters, activity_respiration, basal_respiration, coordinates, irr, k_PAR, surface_PAR, tracers):
    
        # Needs:    Chlorophyll quota (theta_chl in python bfm), Initial slope of PI curve (alpha_chl), Optimal value Epar/EK (p_EpEk_or),
        #           Maximal productivity (rp0), Chl:C relaxation rate (p_tochl_relt)

        # Get carbon concentration
        carbon_index = self.composition.index("c")
        # phyto_carbon = self.conc[carbon_index][iter]
        phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Get chlorophyll concentration
        chl_index = self.composition.index("chl")
        # phyto_chl = self.conc[chl_index][iter]
        phyto_chl = np.array(tracers[self.abbrev].conc[chl_index][iter])

        # Calculate irradiance
        # k_PAR = light_attenuation(parameters, phyto_carbon)
        # irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)

        # x = (self.psn[iter] - self.exu[iter] - activity_respiration) * phyto_carbon
        # y = parameters["initial_PI_slope"] * ( phyto_chl + 1.E-20) * irr
        rho_chl = parameters["chl_quota"] * np.minimum(np.ones_like(self.psn[iter]), (self.psn[iter] - self.exu[iter] - activity_respiration) * phyto_carbon / ( parameters["initial_PI_slope"] * ( phyto_chl + 1.E-20) * irr ))
        chl_opt = parameters["optimal_Epar_Ek"] * parameters["max_photo_rate"] * phyto_carbon / ( parameters["initial_PI_slope"] * irr + 1.E-20 )

        chlorophyll_synthesis = rho_chl * ( self.psn[iter] - self.exu[iter] - activity_respiration ) * phyto_carbon \
                                    - ( self.lys[iter] + basal_respiration ) * phyto_chl - np.maximum(np.zeros_like(phyto_chl), phyto_chl - chl_opt) * parameters["chl_relax_rate"]


        # Update d_dt
        tracers[self.abbrev].d_dt[chl_index] += chlorophyll_synthesis
        # self.d_dt[chl_index] += chlorophyll_synthesis
    

    def exudation(self, iter, parameters, c, p, tracers):
    
        # Locate index of carbon constituent 
        carbon_index_phyto = self.composition.index("c")
        carbon_index_om = tracers[p[0]].composition.index("c")

        # Calculate exudation
        activity = self.psn[iter] * parameters["excreted_fraction"]
        nutrient_stress = self.psn[iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_limitation_factor )
        exudation = (activity + nutrient_stress) * np.array(tracers[self.abbrev].conc[carbon_index_phyto][iter])

        # Update exudation variable
        self.exu[iter] = activity + nutrient_stress

        # Update d_dt
        tracers[c[0]].d_dt[carbon_index_phyto] -= exudation
        tracers[p[0]].d_dt[carbon_index_om] += exudation


    def gross_primary_production(self, iter, parameters, c, p, tracers):
        
        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        # phyto_carbon = self.conc[carbon_index][iter]
        phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Calculate gross primary production
        gross_primary_production = self.psn[iter] * phyto_carbon

        # Update gross primary production variable
        self.gpp[iter] = gross_primary_production

        # Update d_dt
        self.d_dt[carbon_index] += gross_primary_production
        if "o2" in p:   tracers["o2"].d_dt += gross_primary_production / parameters["mw_carbon"]
        if "co2"in c:   tracers["co2"].d_dt -= gross_primary_production
        
    
    def lysis(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
    
        # Needs:    Nutrient limitation, Half sat for stress lysis (h_pnp from python bfm), Activity respiration fraction (d_P0 from python bfm), 
        #           Extra lysis rate (p_seo from python bfm), Half sat for extra lysis (p_sheo from python bfm)

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")
        phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Calculate carbon ratios
        carbon_index = self.composition.index('c')
        if 'n' in self.composition:
            nitrogen_index = self.composition.index('n')
            # nit_carb = self.conc[nitrogen_index] / self.conc[carbon_index]
            nit_carb = np.array(tracers[self.abbrev].conc[nitrogen_index][iter] / tracers[self.abbrev].conc[carbon_index][iter])
        else:
            nit_carb = np.ones_like(phyto_carbon)

        if 'p' in self.composition: 
            phosphorus_index = self.composition.index('p')
            # phos_carb = self.conc[phosphorus_index] / self.conc[carbon_index]
            phos_carb = np.array(tracers[self.abbrev].conc[phosphorus_index][iter] / tracers[self.abbrev].conc[carbon_index][iter])
        else:
            phos_carb = np.ones_like(phyto_carbon)

        # Extract nutrient quotas
        if "no3" in self.nutrients:
            quota_index = self.nutrient_limitation["nutrients"].index("no3")
            min_nitrogen_quota = self.nutrient_limitation["min_quota"][quota_index]
        else:
            min_nitrogen_quota = 0.

        if "po4" in self.nutrients:
            quota_index = self.nutrient_limitation["nutrients"].index("po4")
            min_phosphorus_quota = self.nutrient_limitation["min_quota"][quota_index]
        else: 
            min_phosphorus_quota = 0.
        
        # Calculate fraction of lysis released to dissolved pool
        min_quota = np.minimum(min_nitrogen_quota/(nit_carb + 1.E-20), min_phosphorus_quota/(phos_carb + 1.E-20))
        apportioning_factor = np.minimum(np.ones_like(min_quota), min_quota)

        # Calculate nutrient stress lysis
        nutrient_stress_lysis = ( parameters["max_stress_lysis_rate"] * parameters["half_sat_stress_lysis"] ) / ( self.nutrient_limitation_factor + parameters["half_sat_stress_lysis"] ) \
                                    + ( parameters["extra_lysis_rate"] * phyto_carbon ) / ( phyto_carbon + parameters["half_sat_extra_lysis"] + 1.E-20 )
        
        # Apportion lysis between organic matter pools based on type == particulate or dissolved
        if parameters["om_type"] == "dissolved": lysis = ( 1 - apportioning_factor ) * nutrient_stress_lysis * phyto_carbon
        elif parameters["om_type"] == "particulate": lysis = apportioning_factor * nutrient_stress_lysis * phyto_carbon

        # Calculate concentration ratio
        # concentration_ratio(iter, ic, tracers[c])
        # concentration_ratio(iter, ip, tracers[p])

        # Update lysis variable
        self.lys[iter] = nutrient_stress_lysis

        # Match phyto and organic matter concentration ratios
        ratios = np.zeros(len(tracers[p].conc_ratio))
        for const in self.composition:
            if const in tracers[p].composition:
                index_phyto = self.composition.index(const)
                index_om = tracers[p].composition.index(const)
                ratios[index_om] = tracers[c].conc_ratio[index_phyto]

        # Update d_dt
        tracers[c].d_dt -= ec * tracers[c].conc_ratio * lysis
        tracers[p].d_dt += ep * ratios * lysis
        # tracers[c].d_dt -= ec * tracers[c].conc_ratio * lysis
        # # tracers[p].d_dt += ep * tracers[p].conc_ratio * lysis
        # tracers[p].d_dt += ep * tracers[c].conc_ratio * lysis
    

    def net_primary_production(self, iter, tracers):
        
        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        # phyto_carbon = self.conc[carbon_index][iter]
        phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Calculate losses
        # specific_losses = activity_excretion + activity_respiration + basal_respiration + nutrient_stress_excretion + nutrient_stress_lysis
        specific_losses = self.exu[iter] + self.rsp[iter] + self.lys[iter]

        # Calculate net primary production
        self.npp[iter] = np.maximum( np.zeros_like(phyto_carbon), ( self.psn[iter] - specific_losses ) * phyto_carbon )
    

    def photosynthesis(self, iter, parameters, base_element, coordinates, dz, k_PAR, temperature, surface_PAR, tracers):
        
        # Needs:    Nutrient limitation, Light limitation, Temperature regulation factor, Maximal productivity (Vm) at 10C (rp0 from python bfm, p_sum from fortran bfm)

        # Get concentration of base element in phytoplankton
        base_index = self.composition.index(base_element)
        base_conc = np.array(tracers[self.abbrev].conc[base_index][iter])

        # Get carbon concentration
        if "c" in self.composition:
            carbon_index = self.composition.index("c")
            phyto_carbon = self.conc[carbon_index][iter]
        else:
            phyto_carbon = np.ones_like(base_conc)
        
        # Locate index of chlorophyll constituent if present
        if "chl" in self.composition:   
            chl_index = self.composition.index("chl")
            phyto_chl = self.conc[chl_index][iter]
            pl_pc = phyto_chl / phyto_carbon     # Chl:C ratio (used in light limitation)
        else:
            phyto_chl = base_conc
            pl_pc = 1.

        # Calculate maximal productivity
        if parameters["max_photo_rate"] == "variable":
            Vm = max_growth_rate(parameters, temperature)
        else:
            Vm = parameters["max_photo_rate"]
        
        if iter == 360:
            x=1

        # Calculate light limitation
        if parameters["light_limitation"] in ["monod", "platt", "smith"]:
            # k_PAR = light_attenuation(parameters, phyto_chl)
            irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)
            exp, irr, fI = light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm)
        else:
            fI = parameters["light_limitation"]
            irr = 1.E-20

        # Update photosynthesis variable
        # self.psn[iter] = self.nutrient_limitation_factor * self.temp_regulation_factor * Vm * fI

        fpplim = 1. # this will change once silicate is added
        self.psn[iter] = fpplim * self.temp_regulation_factor * Vm * fI

        return fI, irr
    

    def respiration(self, iter, parameters, c, p, tracers):
        """
        Definition:: Calculates phytoplankton respiration
        """

        # Needs:    Activity and Basal respiration rates (gammaP in python bfm), Activity and Nutrient stress excretion

        # Locate index of carbon constituent
        carbon_index = self.composition.index("c")

        # Get carbon concentration
        # phyto_carbon = self.conc[carbon_index][iter]
        phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Respiration
        activity_respiration = parameters["activity_respiration_frac"] * ( self.psn[iter] - self.exu[iter] )
        basal_respiration = self.temp_regulation_factor * parameters["basal_respiration_rate"]
        total_respiration = activity_respiration + basal_respiration

        respiration = total_respiration * phyto_carbon

        # Update d_dt
        tracers[self.abbrev].d_dt[carbon_index] -= respiration
        if "o2" in c:   tracers["o2"].d_dt -= respiration / parameters["mw_carbon"]
        if "co2" in p:  tracers["co2"].d_dt += respiration        

        # Update respiration variable
        self.rsp[iter] = total_respiration

        return activity_respiration, basal_respiration
    

    def uptake(self, iter, parameters, c, p, ec, ep, ic, tracers):
        """
        Can consume single nutrient (Nitrate/Phosphate) or multiple nutrients (Nitrate AND Ammonium if both available)
        Presence of Ammonium affects uptake of Nitrate
        """

        # Locate carbon index and concentration (if present)
        if "c" in self.composition:
            carbon_index = self.composition.index("c")
            # phyto_carbon = self.conc[carbon_index][iter]
            phyto_carbon = np.array(tracers[p[0]].conc[carbon_index][iter])
        else:
            phyto_carbon = 1.

        # # Get concentration of nutrient in phytoplankton
        # p = p[0]
        # ep = ep[p]
        # nutrient_index = list(ep).index(1.)
        # phyto_nutrient = np.array(tracers[p].conc[nutrient_index][iter])

        # Nutrient uptake can go to dissolved nutrient pool if uptake rate is negative
        # Get concentration of nutrient in phytoplankton
        # Get location of dissolved nutrient pool (if necessary)
        if len(p) > 1:
            phy = list(p).index(self.abbrev)
            ephy = ep[self.abbrev]
            phyto_nutrient_index = list(ephy).index(1.)
            phyto_nutrient = np.array(tracers[self.abbrev].conc[phyto_nutrient_index][iter])

            if phy == 0:    i = 1
            else:           i = 0
            om = p[i]
            eom = ep[om]
            om_nutrient_index = list(eom).index(1.)
            om_nutrient = np.array(tracers[om].conc[om_nutrient_index][iter])
        else:
            phy = p[0]
            ephy = ep[phy]
            phyto_nutrient_index = list(ephy).index(1.)
            phyto_nutrient = np.array(tracers[phy].conc[phyto_nutrient_index][iter])

        # Extract dict
        if len(c) > 1:
            # Get concentrations
            # phyto_nitrogen = self.conc[nit_index][iter]
            no3 = np.array(tracers["no3"].conc[0][iter])
            nh4 = np.array(tracers["nh4"].conc[0][iter])
            
            # Calculate preference for Ammonium uptake
            nh4_preference = parameters["half_sat_nh4_preference"] / ( parameters["half_sat_nh4_preference"] + nh4 + 1.E-20)

            # Calculate maximum nitrogen uptake
            max_uptake_no3 = parameters["specific_affinity"] * no3 * phyto_carbon * nh4_preference
            max_uptake_nh4 = parameters["specific_affinity"] * nh4 * phyto_carbon
            max_uptake_DIN = max_uptake_no3 + max_uptake_nh4

            # Extract nutrient quota
            quota_index = self.nutrient_limitation["nutrients"].index("no3")
            nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]

            # Intracellular missing amount of N
            missing_nit = parameters["max_photo_rate"] * self.temp_regulation_factor * ( parameters["luxury_storage"] * nutrient_quota * phyto_carbon - phyto_nutrient )
            
            # N uptake based on net assimilation of C
            assim_uptake = parameters["luxury_storage"] * nutrient_quota * self.npp[iter]

            # Actual uptake of nitrogen
            actual_uptake = np.minimum(max_uptake_DIN, missing_nit + assim_uptake)

            upt_switch = switch(actual_uptake)

            no3_uptake = upt_switch * actual_uptake * max_uptake_no3 / (max_uptake_DIN + 1.E-20)
            nh4_uptake = upt_switch * actual_uptake * max_uptake_nh4 / (max_uptake_DIN + 1.E-20)

            # phyto_uptake = -actual_uptake * (1. - upt_switch)
            uptake_to_phyto = no3_uptake + nh4_uptake
            uptake_to_om = -actual_uptake * (1. - upt_switch)

            # Update d_dt
            tracers["no3"].d_dt -= no3_uptake
            tracers["nh4"].d_dt -= nh4_uptake
            tracers[self.abbrev].d_dt[phyto_nutrient_index] += uptake_to_phyto
            if len(p) > 1:  
                tracers[self.abbrev].d_dt -= ephy * uptake_to_om
                tracers[om].d_dt += eom * uptake_to_om

            self.uptn[iter] = uptake_to_phyto

        else:
            # Get concentration of nutrient
            c = c[0]
            ec = ec[c]
            ic = ic[c]
            nutrient = np.array(tracers[c].conc[ic][iter])

            # Calculate maximum nutrient uptake
            max_uptake = parameters["specific_affinity"] * nutrient * phyto_carbon

            # Extract nutrient quota
            quota_index = self.nutrient_limitation["nutrients"].index(c)
            nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]

            # Intracellular missing amount of nutrient
            missing = parameters["max_photo_rate"] * self.temp_regulation_factor * ( parameters["luxury_storage"] * nutrient_quota * phyto_carbon - phyto_nutrient )

            # Nutrient uptake based on net assimilation of C
            assim_uptake = parameters["luxury_storage"] * nutrient_quota * self.npp[iter]

            # Actual uptake of nutrient
            actual_uptake = np.minimum(max_uptake, missing + assim_uptake)

            upt_switch = switch(actual_uptake)

            uptake_to_phyto = upt_switch * actual_uptake
            uptake_to_om = -actual_uptake * (1. - upt_switch)

            # Update d_dt
            # tracers[c].d_dt -= np.array(ec) * uptake
            # tracers[self.abbrev].d_dt += ep * phyto_uptake
            tracers[c].d_dt -= np.array(ec) * uptake_to_phyto
            tracers[self.abbrev].d_dt += ephy * uptake_to_phyto
            if len(p) > 1:  
                tracers[self.abbrev].d_dt -= ephy * uptake_to_om
                tracers[om].d_dt += eom * uptake_to_om

            if c == 'po4':  self.uptp[iter] = uptake_to_phyto


    # def uptake(self, iter, parameters, c, p, ec, ep, ic, tracers, surface_PAR):

    #     # light limitation (platt)
    #     fI = parameters["Vm"] * ( 1. - np.exp(-parameters["alpha"] * surface_PAR / ( parameters["Vm"] + 1.E-20 )) )
        
    #     if len(c) > 1:
    #         no3 = np.array(tracers["no3"].conc[0][iter])
    #         nh4 = np.array(tracers["nh4"].conc[0][iter])

    #         f_nh4 = nutrient_limitation(nh4, parameters["kNH4"])
    #         f_no3 = nutrient_limitation(no3, parameters["kNO3"])

    #         uptake_nh4 = parameters["Vm"] * fI * self.temp_regulation_factor * f_nh4
    #         uptake_no3 = parameters["Vm"] * fI * self.temp_regulation_factor * (1. - f_nh4) * f_no3

    #         tracers["no3"].d_dt -= uptake_no3
    #         tracers["nh4"].d_dt -= uptake_nh4
    #         tracers[p].d_dt += uptake_no3 + uptake_nh4

    #     else:
    #         nut = np.array(tracers[c[0]]).conc[0][iter]

    #         f_nut = nutrient_limitation(nut, parameters["kN"])

    #         uptake = parameters["Vm"] * fI * self.temp_regulation_factor * f_nut

    #         tracers[c[0]].d_dt -= uptake
    #         tracers[p[0]].d_dt += ep[p[0]] * uptake

    
    # def photosynthesis(self)





    # def chlorophyll_synthesis(self, iter, parameters, coordinates, surface_PAR):
    #     """
    #     Definition:: Calculates chlorophyll synthesis based on gains and losses of phytoplankton carbon
    #     """
    #     # Locate indexes of carbon and chlorophyll constituents
    #     carbon_index = self.composition.index("c")
    #     chl_index = self.composition.index("chl")

    #     # Get carbon and chlorophyll concentrations
    #     phyto_c = self.conc[carbon_index][iter]
    #     phyto_chl = self.conc[chl_index][iter]

    #     # Calculate chlorophyll-carbon ratio
    #     pl_pc = phyto_chl / phyto_c

    #     # Calculate irradiance
    #     k_PAR = light_attenuation(parameters, phyto_c)
    #     irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)

    #     # Calculate chlorophyll regulation (rho_chl)
    #     factor = ( 1 - parameters["activity_respiration_frac"] ) / ( parameters["initial_PI_slope"] * irrad * phyto_chl + 1E-20 )    # 1E-20 to prevent divide by zero error
    #     rho_chl = parameters["max_chl_c"] * np.minimum(np.ones_like(phyto_chl), np.maximum(np.zeros_like(phyto_chl),factor * (self.gpp[iter] - self.exu[iter])))
    #     # rho_chl = parameters["max_chl_c"] * np.minimum(np.ones_like(phyto_chl), factor * (self.gpp[iter] - self.exu[iter]))

    #     # Calculate chlorophyll synthesis
    #     synthesis = rho_chl * ( 1 - parameters["activity_respiration_frac"] ) * ( self.gpp[iter] - self.exu[iter] ) - pl_pc * ( self.lys[iter] + self.rsp[iter] )
    #     # synthesis = np.maximum(np.zeros_like(phyto_chl),synthesis)  # synthesis >= 0

    #     # Update d_dt
    #     self.d_dt[chl_index] += synthesis


    
    # def exudation(self, iter, parameters, c, p, tracers):
    #     """
    #     Definition:: Calculates the amount of unassimilated carbon released by phytoplankton into dissolved carbon pool
    #     """
    #     # Locate index of carbon constituent 
    #     carbon_index_phyto = self.composition.index("c")
    #     carbon_index_om = tracers[p[0]].composition.index("c")

    #     # Exudation
    #     # exudation = ( parameters["excreted_fraction"] + ( 1 - parameters["excreted_fraction"] ) * ( 1 - self.nutrient_limitation_factor ) ) * self.gpp[iter]
    #     exudation = ( 1 - parameters["excreted_fraction"] ) * ( 1 - self.nutrient_limitation_factor ) * self.gpp[iter]

    #     # Update d_dt
    #     tracers[c[0]].d_dt[carbon_index_phyto] -= exudation
    #     tracers[p[0]].d_dt[carbon_index_om] += exudation

    #     # Update exu variable 
    #     self.exu[iter] = exudation
    

    # def gross_primary_production(self, iter, parameters, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, tracers):
    #     """
    #     Definition:: Calculate gross primary production
    #     """
        
    #     # Get carbon concentration
    #     carbon_index = self.composition.index("c")
    #     phytoc = self.conc[carbon_index][iter]
        
    #     # Locate index of chlorophyll constituent if present
    #     if "chl" in self.composition:   
    #         chl_index = self.composition.index("chl")
    #         phytol = self.conc[chl_index][iter]
    #         pl_pc = phytol / phytoc     # Chl:C ratio (used in light limitation)
    #     else:
    #         pl_pc = 1.

    #     # Calculate growth rate
    #     if parameters["max_photo_rate"] == "variable":
    #         Vm = max_growth_rate(parameters, temperature)
    #     else:
    #         Vm = parameters["max_photo_rate"]

    #     # Calculate light limitation
    #     # if parameters["light_limitation"] == "variable":
    #     if parameters["light_limitation"] in ["monod", "platt", "smith"]:
    #         # k_PAR = light_attenuation(parameters, phytoc)
    #         k_PAR = light_attenuation(parameters, phytol)
    #         irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)
    #         # fI = light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm)
    #         fI = light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm)
    #         # fI = light_limitation(parameters, coordinates, irrad, k_PAR, pl_pc, Vm)
    #     else:
    #         fI = parameters["light_limitation"]

    #     if fI == 1:
    #         x=1
    #     # Calculate gross primary production
    #     gpp = Vm * self.temp_regulation_factor * fI * phytoc

    #     # Update d_dt
    #     self.d_dt[carbon_index] += gpp
    #     if "o2" in p:   tracers["o2"].d_dt += gpp / parameters["mw_carbon"]
    #     if "co2"in c:   tracers["co2"].d_dt -= gpp

    #     # Update gpp variable
    #     self.gpp[iter] = gpp

    #     return irrad


    # def lysis(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
    #     """
    #     Definition:: Calculate lysis rate of phytoplankton to organic matter pools
    #     """
    #     # Extract dict
    #     c = c[0]
    #     p = p[0]
    #     ec = ec[c]
    #     ep = ep[p]
    #     ic = ic[c]
    #     ip = ip[p]
    #     tc = np.array(tracers[c].conc[ic][iter])
    #     tp = np.array(tracers[p].conc[ip][iter])
        
    #     # Calculate carbon ratios
    #     carbon_index = self.composition.index('c')
    #     if 'n' in self.composition:
    #         nitrogen_index = self.composition.index('n')
    #         # nit_carb = self.conc[nitrogen_index] / self.conc[carbon_index]
    #         nit_carb = np.array(tracers[self.abbrev].conc[nitrogen_index][iter] / tracers[self.abbrev].conc[carbon_index][iter])
    #     else:
    #         nit_carb = np.ones_like(tp)

    #     if 'p' in self.composition: 
    #         phosphorus_index = self.composition.index('p')
    #         # phos_carb = self.conc[phosphorus_index] / self.conc[carbon_index]
    #         phos_carb = np.array(tracers[self.abbrev].conc[phosphorus_index][iter] / tracers[self.abbrev].conc[carbon_index][iter])
    #     else:
    #         phos_carb = np.ones_like(tp)
        
    #     # Calculate fraction of lysis released to dissolved pool
    #     min_quota = np.minimum(parameters["min_nitrogen_quota"]/(nit_carb + 1.E-20), parameters["min_phosphorus_quota"]/(phos_carb + 1.E-20))
    #     apportioning_factor = np.minimum(np.ones_like(min_quota), min_quota)

    #     # Calculate nutrient stress limitation
    #     lim = nutrient_limitation(parameters["nutrient_stress_threshold"], self.nutrient_limitation_factor)

    #     # Calculate lysis rate
    #     if parameters["type"] == "dissolved": lysis = ( 1 - apportioning_factor ) * ( lim * parameters["max_lysis_rate"] )
    #     elif parameters["type"] == "particulate": lysis = apportioning_factor * ( lim * parameters["max_lysis_rate"] )

    #     # Calculate concentration ratios
    #     concentration_ratio(iter, ic, tracers[c])
    #     concentration_ratio(iter, ip, tracers[p])

    #     # Update d_dt
    #     tracers[c].d_dt -= ec * tracers[c].conc_ratio * lysis
    #     tracers[p].d_dt += ep * tracers[p].conc_ratio * lysis

    #     # Update lys variable
    #     self.lys += tracers[c].conc_ratio[carbon_index] * lysis


    # def mortality(self, iter, parameters, c, p, ec, ep, ic, ip, tracers):
    #     """
    #     Definition:: Calculates the non-grazing mortality of planktoninc species
    #                  Lysis is parameterized as a quadratic mortality rate (the carbon constituent of lysis is added to the 'lys' variable for the calculation of net primary production)
    #     """
    #     # Extract dict
    #     c = c[0]
    #     p = p[0]
    #     ec = ec[c]
    #     ep = ep[p]
    #     ic = ic[c]
    #     ip = ip[p]
    #     tc = np.array(tracers[c].conc[ic][iter])

    #     # Calculate mortality rate
    #     mortality = ( parameters["mortality_rate"][0] * tc ) + ( parameters["mortality_rate"][1] * (tc**2) )

    #     # Temperature regulation
    #     if self.temp_limited:
    #         mortality = mortality * self.temp_regulation_factor

    #     # Oxygen limitation
    #     if "oxygen_limited" in parameters and parameters["oxygen_limited"]:
    #         fO = nutrient_limitation(tc,parameters["half_sat_oxygen"])
    #         mortality = mortality * (1 - fO)

    #     # Calculate concentration ratios
    #     concentration_ratio(iter, ic, tracers[c])
    #     concentration_ratio(iter, ip, tracers[p])

    #     # Update d_dt
    #     if "partition" in parameters:   # Mortality can be partitioned between dissolved and particulate detrital pools
    #         tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality * parameters["partition"]
    #         tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality * parameters["partition"]

    #     else:
    #         tracers[c].d_dt -= ec * tracers[c].conc_ratio * mortality
    #         tracers[p].d_dt += ep * tracers[p].conc_ratio * mortality

    #     # Update 'lys' variable if phytoplankton contains a carbon constituent
    #     if "c" in self.composition:
    #         carbon_index = self.composition.index("c")
    #         self.lys[iter] += parameters["mortality_rate"][1] * (self.conc[carbon_index][iter]**2)
        

    # def respiration(self, iter, parameters, c, p, tracers):
    #     """
    #     Definition:: Calculates phytoplankton respiration
    #     """
    #     # Locate index of carbon constituent
    #     carbon_index = self.composition.index("c")

    #     # Get carbon concentration
    #     phyto = self.conc[carbon_index][iter]

    #     # Respiration
    #     respiration = parameters["respiration_rate"] * self.temp_regulation_factor * phyto + parameters["activity_respiration_frac"] * (self.gpp[iter] - self.exu[iter])

    #     # Update d_dt
    #     self.d_dt[carbon_index] -= respiration
    #     if "o2" in c:   tracers["o2"].d_dt -= respiration / parameters["mw_carbon"]
    #     if "co2" in p:  tracers["co2"].d_dt += respiration

    #     # Update rsp variable
    #     self.rsp[iter] = respiration


    # def uptake(self, iter, parameters, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, c, p, ec, ep, ic, ip, tracers):
        
    #     # Extract dict
    #     c = c[0]
    #     p = p[0]
    #     ec = ec[c]
    #     ep = ep[p]
    #     ic = ic[c]
    #     ip = ip[p]
    #     tc = np.array(tracers[c].conc[ic][iter])

    #     # Get concentration of nutrient in phytoplankton
    #     nutrient_index = list(ep).index(1.)
    #     tp = np.array(tracers[p].conc[nutrient_index][iter])

    #     # Get carbon concentration
    #     if "c" in self.composition:
    #         carbon_index = self.composition.index("c")
    #         phytoc = self.conc[carbon_index][iter]
    #     else:
    #         phytoc = np.ones_like(tp)
        
    #     # Locate index of chlorophyll constituent if present
    #     if "chl" in self.composition:   
    #         chl_index = self.composition.index("chl")
    #         phyto = self.conc[chl_index][iter]
    #         pl_pc = phyto / phytoc     # Chl:C ratio (used in light limitation)
    #     else:
    #         phyto = tp
    #         pl_pc = 1.

    #     # Calculate growth rate
    #     if parameters["max_photo_rate"] == "variable":
    #         Vm = max_growth_rate(parameters, temperature)
    #     else:
    #         Vm = parameters["max_photo_rate"]
        
    #     # Calculate light limitation
    #     # if parameters["light_limitation"] == "variable":
    #     if parameters["light_limitation"] in ["monod", "platt", "smith"]:
    #         # k_PAR = light_attenuation(parameters, phytoc)
    #         k_PAR = light_attenuation(parameters, phyto)
    #         irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)
    #         # fI = light_limitation(parameters, dz, irrad, k_PAR, mixed_layer_depth, surface_PAR, Vm)
    #         fI = light_limitation(parameters, dz, irrad, k_PAR, pl_pc, Vm)
    #     else:
    #         fI = parameters["light_limitation"]

    #     # Calculate nutrient limitation
    #     if parameters["nutrient_limitation"] == "variable":
    #         fN = nutrient_limitation(tc, parameters["half_sat_nutrient"])
    #     else:
    #         fN = parameters["nutrient_limitation"]
        
    #     # Calculate nutrient uptake
    #     uptake = Vm * self.temp_regulation_factor * fN * fI * tp

    #     # Calculate concentration ratio
    #     # concentration_ratio(iter, ic, tracers[c])
    #     # concentration_ratio(iter, ip, tracers[p])
        
    #     # Update d_dt
    #     tracers[c].d_dt -= np.array(ec) * uptake
    #     tracers[p].d_dt += ep * uptake
    #     # tracers[c].d_dt -= ec * tracers[c].conc_ratio * uptake
    #     # tracers[p].d_dt += ep * tracers[p].conc_ratio * uptake
        
    #     return fI
    

