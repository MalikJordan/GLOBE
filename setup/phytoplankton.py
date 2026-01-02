import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutrient_limitation, monod, temperature_dependence, tracer_elements, switch
from fractions import Fraction
class Phytoplankton():
    """
    
    """

    def __init__(self, abbrev, iters, reactions, **tracer):
        self.abbrev = abbrev
        self.name = tracer["long_name"]
        self.type = tracer["type"]

        # Nutrient limitation
        self.nutrient_limitation = tracer["parameters"]["nutrient_limitation"]
        self.nutrient_limitation_factor = {}
        self.nutrient_colimitation = 0.
        
        # Intracellular nutrinet quotas
        self.cell_quota = tracer["parameters"]["cell_quota"]

        # Light limitation
        if "light_attenuation" in tracer["parameters"]:
            self.light_attenuation = tracer["parameters"]["light_attenuation"]
        else:
            self.light_attenuation = 0.

        # Temperature regulation
        self.temperature_regulation = tracer["parameters"]["temperature_regulation"]
        self.temp_regulation_factor = 1.

        # Oxygen inhibition
        if "oxygen_inhibition" in tracer["parameters"]:
            self.oxygen_inhibition = tracer["parameters"]["oxygen_inhibition"]
            self.oxy_limitation_factor = 1.

        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Phytoplankton: Element required for " + self.name + ". Check documentation and edit input file.")
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
        self.conc_ratio = np.copy(self.cell_quota["opt"])

        # Production 
        # self.exu = np.ones_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.exu = np.zeros_like(self.conc[0,...],dtype=float)   # Exudation (Initialzied to 1 for use in respiration)
        self.gpp = np.zeros_like(self.conc[0,...],dtype=float)   # Gross Primary Production (Initialized to 1 for use in exudation and respiration)
        self.lys = np.zeros_like(self.conc[0,...],dtype=float)  # Lysis (carbon)
        self.npp = np.zeros_like(self.conc[0,...],dtype=float)  # Net Primary Production
        self.psn = np.zeros_like(self.conc[0,...],dtype=float)  # Photosynthesis
        self.rsp = np.zeros_like(self.conc[0,...],dtype=float)  # Respiration
        self.upt = {}   # Uptake

        # self.uptn = np.zeros_like(self.conc[0,...],dtype=float) # Nitrogen uptake
        # self.uptp = np.zeros_like(self.conc[0,...],dtype=float) # Phosophorus uptake

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
        
        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(iter, tracers)

        # Calculate temp regulation factor (if necessary)
        if self.temperature_regulation["temp_limited"]:
            # self.temp_regulation_factor = temperature_dependence(base_temp, temperature, self)
            self.temp_regulation_factor = temperature_dependence(temperature, self)

        check_conc = self.conc[:,iter]

        # Calculate light extiction coefficient
        k_PAR = light_attenuation(self.abbrev, iter, base_element, light_attenuation_water, tracers)
        
        # Calculate rates required for net primary production
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "exudation":                     self.exudation(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "gross_primary_production":      self.gross_primary_production(iter, base_element, reac["parameters"], c, p, tracers)
            if reac["type"] == "lysis":                         self.lysis(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "photosynthesis":                fI, irr = self.photosynthesis(iter, reac["parameters"], coordinates, dz, k_PAR, temperature, surface_PAR)
            if reac["type"] == "respiration":                   activity_respiration, basal_respiration = self.respiration(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)

        # Calculate net primary production
        self.net_primary_production(iter, base_element, tracers)

        # Calculate remaining rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "chlorophyll_synthesis":         
                if self.calc_respiration:   self.chlorophyll_synthesis(iter, base_element, reac["parameters"], activity_respiration, basal_respiration, coordinates, irr, k_PAR, surface_PAR, tracers)
                else:                       self.chlorophyll_synthesis(iter, base_element, reac["parameters"], 0., 0., coordinates, irr, k_PAR, surface_PAR, tracers)
            if reac["type"] == "uptake":    self.uptake(iter, base_element, reac["parameters"], c, p, ec, ep, ic, tracers)


    def add_nutrient(self, nutrient):
        """
        Add nutrients to phytoplankton and append dictionary of uptake rates
        """
        self.upt[nutrient] = np.zeros_like(self.conc[0,...],dtype=float)
        self.nutrient_limitation_factor[nutrient] = np.zeros_like(self.conc[0,...],dtype=float)


    def calculate_nutrient_limitation(self, iter, tracers):
        """
        Definition:: Calculates nutrient limitation factor as either a minimum, product, or sum of all nutrients which limit phytoplankton growth
        """

        fN = []
        for key in self.nutrient_limitation:
            if key != "colimitation" and key != "include":
                # Get nutrient chemical constituent
                element = tracers[key].composition[0]

                # Get index of nutrient chemcical constituent in phytoplankton composition dictionary
                element_index = self.composition.index(element)

                # Get index of nutrient chemcical constituent in cell quota dictionary
                # quota_index = self.cell_quota.index(element)
                quota_index = self.cell_quota["constituents"].index(element)
                
                if self.nutrient_limitation[key]["type"] == "internal":
                    # Calculate nutrient limitation factor
                    func = ( self.conc_ratio[element_index] - self.cell_quota["min"][quota_index]) / ( self.cell_quota["opt"][quota_index] - self.cell_quota["min"][quota_index] )
                    
                    # Ensures nonzero value
                    func = np.maximum(1.E-20*np.ones_like(func), func)

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = func

                    # Append fN for colimitation calculation
                    # if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                    if key in self.nutrient_limitation["include"]:
                        fN.append(func)
                
                elif self.nutrient_limitation[key]["type"] == "external":
                    # Determine if Hill exponent exists for monod function
                    if "exponent" in self.nutrient_limitation[key]:
                        exponent = self.nutrient_limitation[key]["exponent"]
                    else:   # Default to 1. (no scaling)
                        exponent = 1.

                    # Calculate nutrient limitation factor
                    func = monod(tracers[key].conc[...,iter], self.nutrient_limitation[key]["half_sat"], exponent)
                    
                    # Ensures nonzero value
                    func = np.maximum(1.E-20*np.ones_like(func), func)

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = func

                    # Append fN for colimitation calculation
                    # if key in self.nutrient_limitation["colimitation"]["nutrients"]:
                    if "colimitation" in self.nutrient_limitation and key in self.nutrient_limitation["include"]:
                        fN.append(func)

                else:
                    sys.exit("Nutrient limitation type not recognized. Check documentation and edit input file.")

        if "colimitation" in self.nutrient_limitation.keys():   # Multiple nutrients available
            fN = np.array(fN)
            if self.nutrient_limitation["colimitation"] == "minimum":
                self.nutrient_colimitation = np.min(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "product":
                self.nutrient_colimitation = np.prod(fN, axis=0)
            elif self.nutrient_limitation["colimitation"] == "sum":
                self.nutrient_colimitation == np.sum(fN,  axis=0)
            else:
                sys.exit("Nutrient colimitation not recognized. Check documentation and edit input file.")
        else:   # Only one nutrient available
            self.nutrient_colimitation = fN


    def aggregation(self, iter, base_element, parameters, c, ec, ic, tracers):
        """
        Definition:: Calculates the aggregation loss of phytoplankton during high biomass, bloom conditions
        Parameterized as a quadratic (density dependent) loss term
        """
        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Get concentration of base element
        index = self.composition.index(base_element)
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Calculate growth ratio
        growth_ratio = np.minimum(np.ones_like(self.psn[iter]), self.psn[iter]/ ( self.temp_regulation_factor * parameters["aggregation_frac"] * parameters["max_photo_rate"]))

        # Calculate aggregation limit
        agg_limit = (1. - growth_ratio)**2

        # Calculate aggregation loss
        agg_loss = agg_limit * parameters["aggegation_loss"] * (phyto**2)

        # Convert aggregation rate (if necessary)
        if "convert_aggregation" in parameters:
            if parameters["convert_aggregation"] == "cell_quota":
                quota_index = self.nutrient_limitation["nutrients"].index(c)
                agg_loss *= self.nutrient_limitation_factor[quota_index]
            else:
                if isinstance(parameters["convert_aggregation"],(int,float)) and not isinstance(parameters["convert_aggregation"],bool):
                    agg_loss *= parameters["convert_aggregation"]   
                elif isinstance(parameters["convert_aggregation"],str):   
                    agg_loss *= float(Fraction(parameters["convert_aggregation"]))

        # Calculate concentration ratio
        concentration_ratio(iter, index, self)

        # Update d_dt
        self.d_dt -= ec * tracers[c].conc_ratio * agg_loss

        # Apply partition to organic matter group (if necessary)
        if "partition" in parameters:   tracers[p].d_dt += ep * agg_loss * parameters["partition"]
        else:                           tracers[p].d_dt += ep * agg_loss


    def chlorophyll_synthesis(self, iter, base_element, parameters, activity_respiration, basal_respiration, coordinates, irr, k_PAR, surface_PAR, tracers):
    
        # Needs:    Chlorophyll quota (theta_chl in python bfm), Initial slope of PI curve (alpha_chl), Optimal value Epar/EK (p_EpEk_or),
        #           Maximal productivity (rp0), Chl:C relaxation rate (p_tochl_relt)

        # # Get carbon concentration
        # carbon_index = self.composition.index("c")
        # phyto_carbon = np.array(tracers[self.abbrev].conc[carbon_index][iter])

        # Get concentration of base element
        index = self.composition.index(base_element)
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Get chlorophyll concentration
        chl_index = self.composition.index("chl")
        phyto_chl = np.array(tracers[self.abbrev].conc[chl_index][iter])

        # Calculate irradiance
        rho_chl = parameters["chl_quota"] * np.minimum(np.ones_like(self.psn[iter]), (self.psn[iter] - self.exu[iter] - activity_respiration) * phyto / ( parameters["initial_PI_slope"] * ( phyto_chl + 1.E-20) * irr ))
        chl_opt = parameters["optimal_Epar_Ek"] * parameters["max_photo_rate"] * phyto / ( parameters["initial_PI_slope"] * irr + 1.E-20 )

        chlorophyll_synthesis = rho_chl * ( self.psn[iter] - self.exu[iter] - activity_respiration ) * phyto \
                                    - ( self.lys[iter] + basal_respiration ) * phyto_chl - np.maximum(np.zeros_like(phyto_chl), phyto_chl - chl_opt) * parameters["chl_relax_rate"]


        # Update d_dt
        tracers[self.abbrev].d_dt[chl_index] += chlorophyll_synthesis
        # self.d_dt[chl_index] += chlorophyll_synthesis


    def exudation(self, iter, base_element, parameters, c, p, ec, ep, ic, ip, tracers):

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        if parameters["method"] == "constant":
            exudation = parameters["excreted_fraction"] * tracers[c].conc[ic][iter]

        elif parameters["method"] == "photosynthesis":
            # Calculate activity and nutrient stress components
            activity = self.psn[iter] * parameters["excreted_fraction"]
            # nutrient_stress = self.psn[iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_limitation_factor )
            nutrient_stress = self.psn[iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_colimitation )
            
            # Calculate exudation rate
            exudation = (activity + nutrient_stress) * np.array(tracers[self.abbrev].conc[ic][iter])

            # Update exudation variable
            # if ec == base_element:
            if tracers[c].composition[ic] == base_element:
                self.exu[iter] = activity + nutrient_stress

        elif parameters["method"] == "uptake":
            # Calculate total uptake rate for element
            uptake = np.zeros_like(tracers[p].conc[ip][iter])
            for nut in self.upt:
                if tracers[nut].composition[0] == ep:
                    uptake += self.upt[nut][iter]
            
            # Calculate exudation rate
            exudation = parameters["excreted_fraction"] * np.maximum(np.zeros_like(uptake), uptake)

            # Update exudation variable
            if ec == base_element:
                self.exu[iter] = exudation

        # Update d_dt
        # tracers[c[0]].d_dt[ic] -= exudation
        # tracers[p[0]].d_dt[ip] += exudation
        tracers[c].d_dt[ic] -= exudation
        tracers[p].d_dt[ip] += exudation


    def gross_primary_production(self, iter, base_element, parameters, c, p, tracers):
        
        # Locate index of base element
        index = self.composition.index(base_element)

        # Get base element concentration
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Calculate gross primary production
        gross_primary_production = self.psn[iter] * phyto

        # Update gross primary production variable
        self.gpp[iter] = gross_primary_production

        # Update d_dt
        self.d_dt[index] += gross_primary_production
        if "o2" in p:   
            if isinstance(parameters["convert_o2"],(int,float)) and not isinstance(parameters["convert_o2"],bool):
                tracers["o2"].d_dt += gross_primary_production * parameters["convert_o2"]
            elif isinstance(parameters["convert_o2"],str):
                tracers["o2"].d_dt += gross_primary_production * float(Fraction(parameters["convert_o2"]))
        if "co2"in c:   
            if base_element == "c":     tracers["co2"].d_dt -= gross_primary_production
            else:
                if isinstance(parameters["convert_co2"],(int,float)) and not isinstance(parameters["convert_co2"],bool):
                    tracers["co2"].d_dt -= gross_primary_production * parameters["convert_co2"]
                elif isinstance(parameters["convert_co2"],str):   
                    tracers["co2"].d_dt -= gross_primary_production  * float(Fraction(parameters["convert_o2"]))
                

    def lysis(self, iter, base_element, parameters, c, p, ec, ep, ic, ip, tracers):
        
        # Needs:    Nutrient limitation, Half sat for stress lysis (h_pnp from python bfm), Activity respiration fraction (d_P0 from python bfm), 
        #           Extra lysis rate (p_seo from python bfm), Half sat for extra lysis (p_sheo from python bfm)

        # Extract dict
        c = c[0]
        p = p[0]
        ec = ec[c]
        ep = ep[p]
        ic = ic[c]
        ip = ip[p]

        # Locate index of base element
        index = self.composition.index(base_element)
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        if parameters["method"] == "cell_quota":
            # Calculate element ratios
            if 'n' in self.composition:
                nitrogen_index = self.composition.index('n')
                nit_base = np.array(tracers[self.abbrev].conc[nitrogen_index][iter] / tracers[self.abbrev].conc[index][iter])
            else:
                nit_base = np.ones_like(phyto)

            if 'p' in self.composition: 
                phosphorus_index = self.composition.index('p')
                phos_base = np.array(tracers[self.abbrev].conc[phosphorus_index][iter] / tracers[self.abbrev].conc[index][iter])
            else:
                phos_base = np.ones_like(phyto)

            # Extract nutrient quotas
            # if "no3" in self.nutrients:
            #     quota_index = self.nutrient_limitation["nutrients"].index("no3")
            #     min_nitrogen_quota = self.nutrient_limitation["min_quota"][quota_index]
            if "no3" in self.nutrient_limitation:
                quota_index = self.cell_quota["constituents"].index('n')
                min_nitrogen_quota = self.cell_quota["min"][quota_index]
            else:
                min_nitrogen_quota = 0.

            # if "po4" in self.nutrients:
            #     quota_index = self.nutrient_limitation["nutrients"].index("po4")
            #     min_phosphorus_quota = self.nutrient_limitation["min_quota"][quota_index]
            if "po4" in self.nutrient_limitation:
                quota_index = self.cell_quota["constituents"].index('p')
                min_phosphorus_quota = self.cell_quota["min"][quota_index]
            else: 
                min_phosphorus_quota = 0.
            
            # Calculate fraction of lysis released to dissolved pool
            min_quota = np.minimum(min_nitrogen_quota/(nit_base + 1.E-20), min_phosphorus_quota/(phos_base + 1.E-20))
            apportioning_factor = np.minimum(np.ones_like(min_quota), min_quota)

            # Calculate nutrient stress lysis
            # nutrient_stress_lysis = ( parameters["max_stress_lysis_rate"] * parameters["half_sat_stress_lysis"] ) / ( self.nutrient_limitation_factor + parameters["half_sat_stress_lysis"] ) \
            nutrient_stress_lysis = ( parameters["max_stress_lysis_rate"] * parameters["half_sat_stress_lysis"] ) / ( self.nutrient_colimitation + parameters["half_sat_stress_lysis"] ) \
                                        + ( parameters["extra_lysis_rate"] * phyto ) / ( phyto + parameters["half_sat_extra_lysis"] + 1.E-20 )
            
            # Apportion lysis between organic matter pools based on type == particulate or dissolved
            if parameters["om_type"] == "dissolved": lysis = ( 1 - apportioning_factor ) * nutrient_stress_lysis * phyto
            elif parameters["om_type"] == "particulate": lysis = apportioning_factor * nutrient_stress_lysis * phyto

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
        
        elif parameters["method"] == "constant":
            lysis = parameters["lysis_rate"] * self.temp_regulation_factor * ( phyto**2 )

            # Convert lysis rate (if necessary)
            if "convert_lysis" in parameters:
                if parameters["convert_lysis"] == "cell_quota":
                    quota_index = self.nutrient_limitation["nutrients"].index(c)
                    lysis *= self.nutrient_limitation_factor[quota_index]
                else:
                    if isinstance(parameters["convert_lysis"],(int,float)) and not isinstance(parameters["convert_lysis"],bool):
                        lysis *= parameters["convert_lysis"]
                    elif isinstance(parameters["convert_lysis"],str):
                        lysis *= float(Fraction(parameters["convert_lysis"]))
            
            # Update d_dt
            tracers[c].d_dt -= lysis

            # Apply partition to organic matter group (if necessary)
            if "partition" in parameters:   tracers[p].d_dt += lysis * parameters["partition"]
            else:                           tracers[p].d_dt += lysis
    

    def net_primary_production(self, iter, base_element, tracers):
        
        # Locate index of base element
        index = self.composition.index(base_element)

        # Get base element concentration
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Calculate losses
        specific_losses = self.exu[iter] + self.rsp[iter] + self.lys[iter]

        # Calculate net primary production
        self.npp[iter] = np.maximum( np.zeros_like(phyto), ( self.psn[iter] - specific_losses ) * phyto )
    

    def photosynthesis(self, iter, parameters, coordinates, dz, k_PAR, temperature, surface_PAR):
        
        # Maximal productivity
        if parameters["max_photo_rate"] == "eppley":
            Vm = max_growth_rate(parameters, temperature)
        else:
            Vm = parameters["max_photo_rate"]
        
        # Light limitation
        if parameters["light_limitation"] in ["geider", "monod", "platt", "smith"]:
            # Calculate irradiance at surface
            irrad = irradiance(parameters["eps_PAR"], surface_PAR, coordinates, k_PAR)

            # Calculate light limitation
            exp, irr, fI = light_limitation(self, iter, parameters, dz, irrad, k_PAR, Vm)
        else:
            fI = parameters["light_limitation"]
            irr = 1.E-20

        # Update photosynthesis variable
        # self.psn[iter] = self.nutrient_limitation_factor * self.temp_regulation_factor * Vm * fI

        #
        #
        # !!! ADD IN SILICATE LIMITATION OPTION FROM BFM !!!
        #
        #
        fpplim = 1. # this will change once silicate is added
        self.psn[iter] = fpplim * self.temp_regulation_factor * Vm * fI

        return fI, irr
    

    def respiration(self, iter, base_element, parameters, c, p, ec, ep, ic, ip, tracers):
        """
        Definition:: Calculates phytoplankton respiration
        """

        # Needs:    Activity and Basal respiration rates (gammaP in python bfm), Activity and Nutrient stress excretion
        # Extract dict
        if p[0] is not None:
            if len(p) > 1:  # if "co2" also included as produced element
                for index in len(p):
                    if p[index] == "co2":   pass
                    else:   p = p[index]
            else:   p = p[0]
            ep = ep[p]
            ip = ip[p]

        # Locate index of base element
        index = self.composition.index(base_element)

        # Get base_element concentration
        phyto = np.array(tracers[self.abbrev].conc[index][iter])

        # Respiration
        activity_respiration = parameters["activity_respiration_frac"] * ( self.psn[iter] - self.exu[iter] )
        basal_respiration = self.temp_regulation_factor * parameters["basal_respiration_rate"]
        total_respiration = activity_respiration + basal_respiration

        respiration = total_respiration * phyto

        # Update d_dt
        tracers[self.abbrev].d_dt[index] -= respiration
        if p[0] is not None:    tracers[p].d_dt[ip] += respiration
        if "o2" in c:   
            if isinstance(parameters["convert_o2"],(int,float)) and not isinstance(parameters["convert_o2"],bool):
                tracers["o2"].d_dt -= respiration * parameters["convert_o2"]
            elif isinstance(parameters["convert_o2"],str):
                tracers["o2"].d_dt -= respiration * float(Fraction(parameters["convert_o2"]))

        if "co2" in p:  
            if base_element == "c":     tracers["co2"].d_dt += respiration  
            else:                       
                if isinstance(parameters["convert_co2"],(int,float)) and not isinstance(parameters["convert_co2"],bool):
                    tracers["co2"].d_dt += respiration  * parameters["convert_co2"]   
                elif isinstance(parameters["convert_co2"],str):   
                    tracers["co2"].d_dt += respiration  * float(Fraction(parameters["convert_co2"]))

        # Update respiration variable
        self.rsp[iter] = total_respiration

        return activity_respiration, basal_respiration


    def uptake(self, iter, base_element, parameters, c, p, ec, ep, ic, tracers):
        
        # Extract dict
        # ec = ec[c]
        # ep = ep[p]
        # ic = ic[c]
        # ip = ip[p]

        # Identify the chemical constituent of the nutrient(s)
        # element = ec[c[0]]
        element = tracers[c[0]].composition[0]

        # Get concentration of constituent in phytoplankton if present
        if element in self.composition:     element_index = self.composition.index(element)

        if parameters["strategy"] == "coupled":
            coupled_uptake = parameters["coupled_uptake"]
            linked_nutrients = coupled_uptake["link"]

            # Extract uptake rates of linked nutrients
            uptake_rates = []
            for nut in linked_nutrients:
                uptake_rates.append(self.uptake_rates[nut])

            # Calculate total linked uptake rate if multiple linked nutrients are used
            if len(linked_nutrients) > 1:   # use numpy "maximum" to ensure minimum uptake of 0.
                if coupled_uptake["method"] == "max":       linked_uptake = np.maximum(np.maximum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "min":     linked_uptake = np.maximum(np.minimum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "product": linked_uptake = np.maximum(np.prod(uptake_rates), np.zeros_like(self.uptake_rates[nut]))
                elif coupled_uptake["method"] == "sum":     linked_uptake = np.maximum(np.sum(uptake_rates), np.zeros_like(self.uptake_rates[nut]))

            if isinstance(parameters["convert_uptake"],(int,float)) and not isinstance(parameters["convert_uptake"],bool):
                uptake = self.nutrient_limitation_factor[element] * linked_uptake * coupled_uptake["convert_uptake"]
            elif isinstance(parameters["convert_uptake"],str):
                uptake = self.nutrient_limitation_factor[element] * linked_uptake * float(Fraction(parameters["convert_uptake"]))

            

            # Update d_dt
            tracers[c].d_dt -= np.array(ec) * np.maximum(uptake, np.zeros_like(uptake))
            if element in self.composition:     self.d_dt[element_index] += np.maximum(uptake, np.zeros_like(uptake))

        elif parameters["strategy"] == "independent":

            if len(c) > 1:  # Used if no3 and nh4 are consumed together rather than individually
                if len(p) > 1:  # If uptake can be source of organic matter
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

                # Get concentration of element in phytoplankton
                if parameters["form"] == "affinity":    # use specific affinity
                    index = self.composition.index(base_element)
                else:   # use nutrient constituent
                    index = self.composition.index("n")
                # if "n" in self.composition:
                #     index = self.composition.index("n")
                # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                #     index = self.composition.index(base_element)
                phyto = np.array(self.conc[index][iter])

                # Get nutrient concentrations
                no3 = np.array(tracers["no3"].conc[0][iter])
                nh4 = np.array(tracers["nh4"].conc[0][iter])
                
                # Calculate preference for Ammonium uptake
                nh4_preference = parameters["half_sat_nh4_preference"] / ( parameters["half_sat_nh4_preference"] + nh4 + 1.E-20)

                # Calculate maximum nitrogen uptake
                max_uptake_no3 = parameters["specific_affinity"] * no3 * phyto * nh4_preference
                max_uptake_nh4 = parameters["specific_affinity"] * nh4 * phyto
                max_uptake_DIN = max_uptake_no3 + max_uptake_nh4

                # Extract nutrient quota
                # quota_index = self.nutrient_limitation["nutrients"].index("no3")
                # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
                quota_index = self.cell_quota["constituents"].index('n')
                nutrient_quota = self.cell_quota["opt"][quota_index]

                # Intracellular missing amount of N
                missing_nit = parameters["max_photo_rate"] * self.temp_regulation_factor * ( parameters["luxury_storage"] * nutrient_quota * phyto - phyto_nutrient )
                
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

                # Calculate n2 uptake for nitrogen fixers (if necessary)
                if "n2" in c:
                    n2_uptake = ( 1. - self.nutrient_limitation_factor["no3"] - self.nutrient_limitation_factor["nh4"] ) * self.psn[iter] * phyto_nutrient
                    uptake_to_phyto += np.maximum(n2_uptake, np.zeros_like(n2_uptake))

                # Update d_dt
                tracers["no3"].d_dt -= no3_uptake
                tracers["nh4"].d_dt -= nh4_uptake
                tracers[self.abbrev].d_dt[phyto_nutrient_index] += uptake_to_phyto
                if len(p) > 1:  
                    tracers[self.abbrev].d_dt -= ephy * uptake_to_om
                    tracers[om].d_dt += eom * uptake_to_om

                # self.uptn[iter] = uptake_to_phyto

            else:
                c = c[0]
                if c == "no3":
                    # Get concentration of element in phytoplankton
                    if "n" in self.composition:
                        index = self.composition.index("n")
                    else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                        index = self.composition.index(base_element)

                    phyto = np.array(self.conc[index][iter])

                    # Determine uptake strategy
                    if parameters["basis"] == "constant":
                        uptake = parameters["constant"] * self.nutrient_limitation_factor["no3"] * phyto
                    elif parameters["basis"] == "growth":
                        uptake = self.psn[iter] * self.nutrient_limitation_factor["no3"] * phyto
                    elif parameters["basis"] == "nutrient":
                        pass
                    
                    # Multiply by Monod function of nutrient limitation (if necessary)
                    if "nh4" in tracers and self.nutrient_limitation["no3"]["nh4_inhibited"]:    # no3_lim / (no3_lim + nh4_lim)
                        uptake *= monod(self.nutrient_limitation_factor["no3"], self.nutrient_limitation_factor["nh4"], 1.)

                    # Update d_dt
                    tracers[c].d_dt -= np.maximum(uptake, np.zeros_like(uptake))
                    if "n" in self.composition:     self.d_dt[index] += np.maximum(uptake, np.zeros_like(uptake))

                elif c == "nh4":
                    # Get concentration of element in phytoplankton
                    if "n" in self.composition:
                        index = self.composition.index("n")
                    else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                        index = self.composition.index(base_element)

                    phyto = np.array(self.conc[index][iter])

                    # Determine uptake strategy
                    if parameters["basis"] == "constant":
                        uptake = parameters["constant"] * self.nutrient_limitation_factor["nh4"] * phyto
                    elif parameters["basis"] == "growth":
                        uptake = self.psn[iter] * self.nutrient_limitation_factor["nh4"] * phyto
                    elif parameters["basis"] == "nutrient":
                        pass

                    # Multiply by Monod function of nutrient limitation (if necessary)
                    if self.nutrient_limitation["no3"]["nh4_inhibited"]:    # nh4_lim / (nh4_lim + no3_lim)
                        uptake *= monod(self.nutrient_limitation_factor["nh4"], self.nutrient_limitation_factor["no3"], 1.)

                    # Update d_dt
                    tracers[c].d_dt -= np.maximum(uptake, np.zeros_like(uptake))
                    if "n" in self.composition:     self.d_dt[index] += np.maximum(uptake, np.zeros_like(uptake))

                elif c == "po4":
                    # Get concentration of element in phytoplankton
                    # if "p" in self.composition:
                    #     index = self.composition.index("p")
                    # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                    #     index = self.composition.index(base_element)

                    if parameters["form"] == "affinity": # use specific affinity
                        index = self.composition.index(base_element)
                    else:   # use nutrient constituent
                        index = self.composition.index("p")
                    phyto = np.array(self.conc[index][iter])

                    # Determine uptake strategy
                    if parameters["basis"] == "constant":
                        uptake = parameters["constant"] * self.nutrient_limitation["po4"] * phyto
                    elif parameters["basis"] == "growth":
                        pass
                    elif parameters["basis"] == "nutrient":
                        if len(p) > 1: # If uptake can be source of organic matter
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
                        
                        # Get concentration of nutrient
                        # c = c[0]
                        ec = ec[c]
                        ic = ic[c]
                        nutrient = np.array(tracers[c].conc[ic][iter])

                        # Calculate maximum nutrient uptake
                        max_uptake = parameters["specific_affinity"] * nutrient * phyto

                        # Extract nutrient quota
                        # quota_index = self.nutrient_limitation["nutrients"].index(c)
                        # nutrient_quota = self.nutrient_limitation["opt_quota"][quota_index]
                        quota_index = self.cell_quota["constituents"].index('p')
                        nutrient_quota = self.cell_quota["opt"][quota_index]

                        # Intracellular missing amount of nutrient
                        missing = parameters["max_photo_rate"] * self.temp_regulation_factor * ( parameters["luxury_storage"] * nutrient_quota * phyto - phyto_nutrient )

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

                        # if c == 'po4':  self.uptp[iter] = uptake_to_phyto
                    
                elif c == "fe":
                    # Get concentration of element in phytoplankton
                    if "fe" in self.composition:
                        index = self.composition.index("fe")
                    else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                        index = self.composition.index(base_element)

                    phyto = np.array(self.conc[index][iter])

                    # Determine uptake strategy
                    if parameters["basis"] == "constant":
                        uptake = parameters["constant"] * self.temp_regulation_factor * self.nutrient_limitation["fe"] * phyto
                    elif parameters["basis"] == "growth":
                        pass
                    elif parameters["basis"] == "nutrient":
                        pass
                    
                    # Mulitply by conversion factor (if needed)
                    if "convert_uptake" in parameters:
                        if isinstance(parameters["convert_uptake"],(int,float)) and not isinstance(parameters["convert_uptake"],bool):
                            uptake *= coupled_uptake["convert_uptake"]
                        elif isinstance(parameters["convert_uptake"],str):
                            uptake *= float(Fraction(parameters["convert_uptake"]))

                    # Uptake is zero if Fe:Base ratio meets or exceeds maximum ratio
                    if "fe" in self.composition:    # Only need to calculate if iron is directly resolved
                        concentration_ratio(iter, index, self.conc)
                        if self.conc_ratio[index] >= self.cell_quota["max"]:    uptake = np.zeros_like(uptake)

                    # Update d_dt
                    tracers[c].d_dt -= uptake
                    if "fe" in self.composition:    self.d_dt[index] += uptake

                elif c == "sio4":    
                    pass
        