import os
import sys
import numpy as np
from functions.seasonal_cycling import *
from functions.other_functions import concentration_ratio, irradiance, light_attenuation, light_limitation, max_growth_rate, nutrient_limitation, monod, temperature_dependence, tracer_elements, switch
from fractions import Fraction
from pom.check_phy import nitr_lim, phos_lim, mult_lim, photo
class Phytoplankton():
    """
    
    """

    def __init__(self, abbrev, base_element, iters, num_layers, reactions, **tracer):
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

        # Sedimentation
        if "sedimentation" in tracer["parameters"]:
            if tracer["parameters"]["sedimentation"]["sinking"] == True:
                self.sinking_velocity = np.ones(num_layers-1) * tracer["parameters"]["sedimentation"]["sinking_rate"]
                self.sinking_velocity[-1] = tracer["parameters"]["sedimentation"]["burial_velocity"]

        # Composition and concentration arrays
        self.composition = []
        conc = []
        if len(tracer["composition"]) < 1:
            sys.exit("Phytoplankton: Element required for " + self.name + ". Check documentation and edit input file.")
        else:
            for key in tracer["composition"]:
                available_elements = ['c','n','p','chl','fe','si','caco3']
                if key in available_elements:
                    # Add constituent to composition/concentration
                    self.composition.append(key)

                    # Set initial conditions
                    if isinstance(tracer["composition"][key], str): # Read initial conditions from file
                        conc.append( np.fromfile(os.getcwd() + tracer["composition"][key]) )
                    elif isinstance(tracer["composition"][key], (int,float)): # Create array of initial conditions
                        if num_layers == 1: # 0d configuration
                            conc.append(tracer["composition"][key])
                        else: # 1d configuration
                            conc.append(tracer["composition"][key] * np.ones(num_layers))
                    elif isinstance(tracer["composition"][key], (list,np.ndarray)):
                        conc.append(np.array(tracer["composition"][key]))
                    elif isinstance(tracer["composition"][key], dict): # Create array of initial conditions based off ratio to base element
                        index = self.composition.index(base_element)
                        if num_layers == 1: # 0d configuration
                            conc.append(tracer["composition"][key][base_element] * conc[index] )
                        else: # 1d configuration
                            conc.append(tracer["composition"][key][base_element] * conc[index] * np.ones(num_layers))

                else:
                    sys.exit("Phytoplankton: Element '" + key + "' not recognized. Check documentation and edit input file.")
                
        # hold = np.zeros((len(conc),iters),dtype=np.ndarray)
        # hold[...,0] = conc
        # self.conc = hold
        # self.d_dt = np.zeros_like(conc)
        # self.conc_ratio = np.copy(self.cell_quota["opt"])

        if num_layers > 1:  # Model as "boxes" between layers (num_layers-1)
            self.conc = np.zeros((len(self.composition),num_layers-1,iters),dtype=float)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const][:-1]
        else:   # Model as single box
            self.conc = np.zeros((len(self.composition),iters),dtype=float)
            for const in range(0,len(self.composition)):
                self.conc[const,:,0] = conc[const]
        self.d_dt = np.zeros_like(self.conc[...,0],dtype=float)
        self.conc_ratio = np.ones_like(self.conc[...,0],dtype=float)
        

        # Production 
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


    # def phyto(self, iter, base_element, base_temp, light_attenuation_water, coordinates, dz, mixed_layer_depth, surface_PAR, temperature, tracers):
    def phyto(self, iter, base_element, physical, tracers):
        
        # Calculate nutrient limitation
        self.calculate_nutrient_limitation(iter, tracers)

        if iter < 10:
            dif_nitr_lim = self.nutrient_limitation_factor["no3"] - nitr_lim[iter]
            dif_phos_lim = self.nutrient_limitation_factor["po4"] - phos_lim[iter]
            dif_mult_lim = self.nutrient_colimitation - mult_lim[iter]

        # Calculate temp regulation factor (if necessary)
        if self.temperature_regulation["temp_limited"]:
            self.temp_regulation_factor = temperature_dependence(physical["bgc_phys_vars"]["temperature"], self)

        # check_conc = self.conc[:,iter]

        # Calculate light extiction coefficient
        k_PAR = light_attenuation(self.abbrev, iter, base_element, physical["environment"]["light_attenuation_water"], tracers)
        
        # Calculate rates required for net primary production
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "exudation":                     self.exudation(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            if reac["type"] == "gross_primary_production":      self.gross_primary_production(iter, base_element, reac["parameters"], c, p, tracers)
            if reac["type"] == "lysis":                         self.lysis(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)
            # if reac["type"] == "photosynthesis":                fI, irr = self.photosynthesis(iter, reac["parameters"], physical["bgc_phys_vars"]["z"], physical["bgc_phys_vars"]["dz"], k_PAR, physical["bgc_phys_vars"]["temperature"], physical["bgc_phys_vars"]["surface_PAR"])
            if reac["type"] == "photosynthesis":                fI, irr = self.photosynthesis(iter, reac["parameters"], physical["bgc_phys_vars"]["z"], physical["bgc_phys_vars"]["z"], k_PAR, physical["bgc_phys_vars"]["temperature"], physical["bgc_phys_vars"]["surface_PAR"])
            if reac["type"] == "respiration":                   activity_respiration, basal_respiration = self.respiration(iter, base_element, reac["parameters"], c, p, ec, ep, ic, ip, tracers)

            x=1
        # Calculate net primary production
        self.net_primary_production(iter, base_element, tracers)

        # Calculate remaining rates
        for reac in self.reactions:
            c, p, ec, ep, ic, ip = tracer_elements(base_element, reac, tracers)
            if reac["type"] == "chlorophyll_synthesis":         
                if self.calc_respiration:   self.chlorophyll_synthesis(iter, base_element, reac["parameters"], activity_respiration, basal_respiration, physical["bgc_phys_vars"]["z"], irr, k_PAR, physical["bgc_phys_vars"]["surface_PAR"], tracers)
                else:                       self.chlorophyll_synthesis(iter, base_element, reac["parameters"], 0., 0., physical["bgc_phys_vars"]["z"], irr, k_PAR, physical["bgc_phys_vars"]["surface_PAR"], tracers)
            if reac["type"] == "uptake":    self.uptake(iter, base_element, reac["parameters"], c, p, ec, ep, ic, tracers)
            x=1
        x=1

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
                    # func = np.maximum(1.E-20*np.ones_like(func), func)
                    func = np.minimum(np.ones_like(func),np.maximum(1.E-20*np.ones_like(func), func))   # maximum value of 1.

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = np.minimum(np.ones_like(func), func)

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
                    # func = np.maximum(1.E-20*np.ones_like(func), func)
                    func = np.minimum(np.ones_like(func),np.maximum(1.E-20*np.ones_like(func), func))   # maximum value of 1.

                    # Update dictionary
                    self.nutrient_limitation_factor[key] = np.minimum(np.ones_like(func), func)

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
        phyto = np.array(tracers[self.abbrev].conc[index,:,iter])

        # Get chlorophyll concentration
        chl_index = self.composition.index("chl")
        phyto_chl = np.array(tracers[self.abbrev].conc[chl_index,:,iter])

        # Calculate irradiance
        rho_chl = parameters["chl_quota"] * np.minimum(np.ones_like(self.psn[:,iter]), (self.psn[:,iter] - self.exu[:,iter] - activity_respiration) * phyto / ( parameters["initial_PI_slope"] * ( phyto_chl + 1.E-20) * irr ))
        chl_opt = parameters["optimal_Epar_Ek"] * parameters["max_photo_rate"] * phyto / ( parameters["initial_PI_slope"] * irr + 1.E-20 )

        chlorophyll_synthesis = rho_chl * ( self.psn[:,iter] - self.exu[:,iter] - activity_respiration ) * phyto \
                                    - ( self.lys[:,iter] + basal_respiration ) * phyto_chl - np.maximum(np.zeros_like(phyto_chl), phyto_chl - chl_opt) * parameters["chl_relax_rate"]


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
            exudation = parameters["excreted_fraction"] * tracers[c].conc[ic,:,iter]

        elif parameters["method"] == "photosynthesis":
            # Calculate activity and nutrient stress components
            activity = self.psn[:,iter] * parameters["excreted_fraction"]
            # nutrient_stress = self.psn[iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_limitation_factor )
            nutrient_stress = self.psn[:,iter] * ( 1. - parameters["excreted_fraction"] ) * ( 1. - self.nutrient_colimitation )
            
            # Calculate exudation rate
            exudation = (activity + nutrient_stress) * np.array(tracers[self.abbrev].conc[ic,:,iter])

            # Update exudation variable
            # if ec == base_element:
            if tracers[c].composition[ic] == base_element:
                self.exu[:,iter] = activity + nutrient_stress

        elif parameters["method"] == "uptake":
            # Calculate total uptake rate for element
            uptake = np.zeros_like(tracers[p].conc[ip,:,iter])
            for nut in self.upt:
                if tracers[nut].composition[0] == ep:
                    uptake += self.upt[nut,:,iter]
            
            # Calculate exudation rate
            exudation = parameters["excreted_fraction"] * np.maximum(np.zeros_like(uptake), uptake)

            # Update exudation variable
            if ec == base_element:
                self.exu[:,iter] = exudation

        # Update d_dt
        # tracers[c[0]].d_dt[ic] -= exudation
        # tracers[p[0]].d_dt[ip] += exudation
        tracers[c].d_dt[ic] -= exudation
        tracers[p].d_dt[ip] += exudation


    def gross_primary_production(self, iter, base_element, parameters, c, p, tracers):
        
        # Locate index of base element
        index = self.composition.index(base_element)

        # Get base element concentration
        phyto = np.array(tracers[self.abbrev].conc[index,:,iter])

        # Calculate gross primary production
        gross_primary_production = self.psn[:,iter] * phyto

        # Update gross primary production variable
        self.gpp[:,iter] = gross_primary_production

        # Update d_dt
        self.d_dt[index,:] += gross_primary_production
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
        phyto = np.array(tracers[self.abbrev].conc[index,:,iter])

        if parameters["method"] == "cell_quota":
            # Calculate element ratios
            if 'n' in self.composition:
                nitrogen_index = self.composition.index('n')
                nit_base = np.array(tracers[self.abbrev].conc[nitrogen_index,:,iter] / tracers[self.abbrev].conc[index,:,iter])
            else:
                nit_base = np.ones_like(phyto)

            if 'p' in self.composition: 
                phosphorus_index = self.composition.index('p')
                phos_base = np.array(tracers[self.abbrev].conc[phosphorus_index,:,iter] / tracers[self.abbrev].conc[index,:,iter])
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
            # if parameters["om_type"] == "dissolved": lysis = ( 1 - apportioning_factor ) * nutrient_stress_lysis * phyto
            # elif parameters["om_type"] == "particulate": lysis = apportioning_factor * nutrient_stress_lysis * phyto

            if tracers[p].form == "dissolved": lysis = ( 1 - apportioning_factor ) * nutrient_stress_lysis * phyto
            elif tracers[p].form == "particulate": lysis = apportioning_factor * nutrient_stress_lysis * phyto

            # Update lysis variable
            self.lys[:,iter] = nutrient_stress_lysis

            # Match phyto and organic matter concentration ratios
            # ratios = np.zeros(len(tracers[p].conc_ratio))
            ratios = np.zeros_like(tracers[p].conc_ratio)
            for const in self.composition:
                if const in tracers[p].composition:
                    index_phyto = self.composition.index(const)
                    index_om = tracers[p].composition.index(const)
                    ratios[index_om] = tracers[c].conc_ratio[index_phyto]

            # Update d_dt
            # tracers[c].d_dt -= ec * tracers[c].conc_ratio * lysis
            # tracers[p].d_dt += ep * ratios * lysis
            for i in range(len(ec)):
                tracers[c].d_dt[i] -= ec[i] * tracers[c].conc_ratio[i] * lysis
            for j in range(len(ep)):
                tracers[p].d_dt[j] += ep[j] * ratios[j] * lysis
            # tracers[c].d_dt -= ((tracers[c].conc_ratio).T * ec).T * lysis
            # tracers[p].d_dt += (ratios.T * ep).T * lysis

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
        phyto = np.array(tracers[self.abbrev].conc[index,:,iter])

        # Calculate losses
        specific_losses = self.exu[:,iter] + self.rsp[:,iter] + self.lys[:,iter]

        # Calculate net primary production
        self.npp[:,iter] = np.maximum( np.zeros_like(phyto), ( self.psn[:,iter] - specific_losses ) * phyto )
    

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
        self.psn[:,iter] = fpplim * self.temp_regulation_factor * Vm * fI

        diff_psn = self.psn[0,iter] - photo[iter]
        if abs(diff_psn) > 1E-18:
            x=1

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
        phyto = np.array(tracers[self.abbrev].conc[index,:,iter])

        # Respiration
        activity_respiration = parameters["activity_respiration_frac"] * ( self.psn[:,iter] - self.exu[:,iter] )
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
        self.rsp[:,iter] = total_respiration

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
                    phyto_nutrient = np.array(tracers[self.abbrev].conc[phyto_nutrient_index,:,iter])

                    if phy == 0:    i = 1
                    else:           i = 0
                    om = p[i]
                    eom = ep[om]
                    om_nutrient_index = list(eom).index(1.)
                    om_nutrient = np.array(tracers[om].conc[om_nutrient_index,:,iter])
                else:
                    phy = p[0]
                    ephy = ep[phy]
                    phyto_nutrient_index = list(ephy).index(1.)
                    phyto_nutrient = np.array(tracers[phy].conc[phyto_nutrient_index,:,iter])

                # Get concentration of element in phytoplankton
                if parameters["form"] == "affinity":    # use specific affinity
                    index = self.composition.index(base_element)
                else:   # use nutrient constituent
                    index = self.composition.index("n")
                # if "n" in self.composition:
                #     index = self.composition.index("n")
                # else:   # If element doesn't isn't directly resolved, use base element with conversion factor
                #     index = self.composition.index(base_element)
                phyto = np.array(self.conc[index,:,iter])

                # Get nutrient concentrations
                no3 = np.array(tracers["no3"].conc[0,:,iter])
                nh4 = np.array(tracers["nh4"].conc[0,:,iter])
                
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
                assim_uptake = parameters["luxury_storage"] * nutrient_quota * self.npp[:,iter]

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
                    for i in range(len(ephy)):
                        tracers[self.abbrev].d_dt[i] -= ephy[i] * uptake_to_om
                    for j in range(len(eom)):
                        tracers[om].d_dt[j] += eom[j] * uptake_to_om

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
                    phyto = np.array(self.conc[index,:,iter])

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
                            phyto_nutrient = np.array(tracers[self.abbrev].conc[phyto_nutrient_index,:,iter])

                            if phy == 0:    i = 1
                            else:           i = 0
                            om = p[i]
                            eom = ep[om]
                            om_nutrient_index = list(eom).index(1.)
                            om_nutrient = np.array(tracers[om].conc[om_nutrient_index,:,iter])
                        else:
                            phy = p[0]
                            ephy = ep[phy]
                            phyto_nutrient_index = list(ephy).index(1.)
                            phyto_nutrient = np.array(tracers[phy].conc[phyto_nutrient_index,:,iter])
                        
                        # Get concentration of nutrient
                        # c = c[0]
                        ec = ec[c]
                        ic = ic[c]
                        nutrient = np.array(tracers[c].conc[ic,:,iter])

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
                        assim_uptake = parameters["luxury_storage"] * nutrient_quota * self.npp[:,iter]

                        # Actual uptake of nutrient
                        actual_uptake = np.minimum(max_uptake, missing + assim_uptake)
                        # actual_uptake = np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))

                        upt_switch = switch(actual_uptake)

                        uptake_to_phyto = upt_switch * actual_uptake
                        uptake_to_om = -actual_uptake * (1. - upt_switch)

                        # Update d_dt
                        # tracers[c].d_dt -= np.array(ec) * uptake
                        # tracers[self.abbrev].d_dt += ep * phyto_uptake
                        tracers[c].d_dt -= np.array(ec) * uptake_to_phyto
                        # tracers[self.abbrev].d_dt += ephy * uptake_to_phyto
                        if len(p) > 1:  
                            # tracers[self.abbrev].d_dt -= ephy * uptake_to_om
                            # tracers[om].d_dt += eom * uptake_to_om
                            for i in range(len(ephy)):
                                # tracers[self.abbrev].d_dt[i] += ephy[i] * uptake_to_phyto
                                tracers[self.abbrev].d_dt[i] += ephy[i] *  np.minimum(max_uptake, missing + assim_uptake)
                                # tracers[self.abbrev].d_dt[i] += ephy[i] * np.maximum(np.zeros_like(max_uptake), np.minimum(max_uptake,missing + assim_uptake))
                            for j in range(len(eom)):
                                tracers[om].d_dt[j] += eom[j] * uptake_to_om
                        else:
                            for i in range(len(ephy)):
                                tracers[self.abbrev].d_dt[i] += ephy[i] * uptake_to_phyto
                        

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
        

assim_uptake_po4_fortran = [1.6814743045847532E-004, 1.6931080142746542E-004, 1.7034718892522456E-004, 1.7124839238892624E-004, 1.7200641914569083E-004, 1.7261361994317368E-004, 1.7306281544991105E-004, 1.7334742814187104E-004, 1.7346160323831496E-004, 1.7340031546787010E-004, 1.7315946010880919E-004, 1.7273593534617560E-004, 1.7212769082820340E-004, 1.7133377211496077E-004, 1.7035434040819778E-004, 1.6919067063121392E-004, 1.6784513011071069E-004, 1.6632114248375522E-004, 1.6462313935402135E-004, 1.6275650204654449E-004, 1.6072749261816771E-004, 1.5854317708766078E-004, 1.5621134427765224E-004, 1.5374041409832085E-004, 1.5113933855591651E-004, 1.4841750169954302E-004, 1.4558462443464392E-004, 1.4265067271827456E-004, 1.3962576863043333E-004, 1.3652010582862866E-004, 1.3334386988191963E-004, 1.3010716451778638E-004, 1.2681994397299381E-004, 1.2349195209338499E-004, 1.2013266794871410E-004, 1.1675125841887996E-004, 1.1335653731508610E-004, 1.0995693126504217E-004, 1.0656045181663018E-004, 1.0317467307159795E-004, 9.9806713758458234E-005, 9.6463225389579004E-005, 9.3150386448991637E-005, 8.9873903090643864E-005, 8.6639008396373639E-005, 8.3450467622888839E-005, 8.0312587670042159E-005, 7.7229229436542887E-005, 7.4203822420670947E-005, 7.1239381050164309E-005, 6.8338522726765148E-005, 6.5503486835033779E-005, 6.2736154754105639E-005, 6.0038070212806188E-005, 5.7410460073713753E-005, 5.4854254981585908E-005, 5.2370110003010785E-005, 4.9958424784680614E-005, 4.7619363282449247E-005, 4.5352872445592850E-005, 4.3158699298633454E-005, 4.1036407654583617E-005, 3.8985393924133195E-005, 3.7004904147754821E-005, 3.5094049396605343E-005, 3.3251820012956318E-005, 3.1477099567931615E-005, 2.9768678658235373E-005, 2.8125268217343930E-005, 2.6545512156783862E-005, 2.5027998953392441E-005, 2.3571271936914568E-005, 2.2173838038410357E-005, 2.0834175805129546E-005, 1.9550742750530925E-005, 1.8321982088222062E-005, 1.7146328852375955E-005, 1.6022215500419999E-005, 1.4948076960008753E-005, 1.3922355033718204E-005, 1.2943501438227491E-005, 1.2009981926658952E-005, 1.1120280712209818E-005, 1.0272904163103515E-005, 9.4663837372109869E-006, 8.6992783507400361E-006, 7.9701765038140244E-006, 7.2776985410142118E-006, 6.6204992006639646E-006, 5.9972703866830856E-006, 5.4067439563747782E-006, 4.8476949803307107E-006, 4.3189435741996670E-006, 3.8193575591245758E-006, 3.3478549361397293E-006, 2.9034055154335121E-006, 2.4850307436573477E-006, 2.0918005771205208E-006, 1.7228260615431317E-006, 1.3772473094822884E-006, 1.0520394109618397E-006, 7.4841771380025374E-007, 4.6804596173191955E-007, 2.0961681546458945E-007, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000]
max_uptake_po4_fortran = [2.3322032227480109E-006, 2.2159108475233858E-006, 2.0998711028945595E-006, 1.9842753398030213E-006, 1.8693113766360872E-006, 1.7551634056100194E-006, 1.6420118394842449E-006, 1.5300333072479151E-006, 1.4194006239194659E-006, 1.3102827304881900E-006, 1.2028446317255370E-006, 1.0972473799608120E-006, 9.9364776330898771E-007, 8.9219762015618414E-007, 7.9304387361634361E-007, 6.9632883884215754E-007, 6.0219048337636842E-007, 5.1076244145566155E-007, 4.2217397946064903E-007, 3.3654996770127697E-007, 2.5401084762582901E-007, 1.7467246886914725E-007, 9.8645668974234488E-008, 2.6036578892845950E-008, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 2.5010916228694582E-008, 1.0060151774625405E-007, 1.8012627835063843E-007, 2.6350673229612082E-007, 3.5065963302137560E-007, 4.4149695547807778E-007, 5.3592580033760004E-007, 6.3384839444999480E-007, 7.3516211394711313E-007, 8.3975950362100850E-007, 9.4752830529255471E-007, 1.0583515079658796E-006, 1.1721074286052046E-006, 1.2886698320185539E-006, 1.4079080881666861E-006, 1.5296872801970861E-006, 1.6538682266194152E-006, 1.7803070497872489E-006, 1.9088542237249201E-006, 2.0393541288224435E-006, 2.1716448100857106E-006, 2.3055581376139173E-006, 2.4409212041729019E-006, 2.5775565360750501E-006, 2.7152820982625041E-006, 2.8539115018167596E-006, 2.9932540151126812E-006, 3.1331135258763517E-006, 3.2732874696606639E-006, 3.4135668252029530E-006, 3.5537369395992472E-006, 3.6935801907979173E-006, 3.8328792119465435E-006, 3.9714204567534590E-006, 4.1089993985000854E-006, 4.2454265102210280E-006, 4.3805354554980369E-006, 4.5141929016444797E-006, 4.6463082619130499E-006, 4.7768394966921376E-006, 4.9057870009249371E-006, 5.0331612569504544E-006, 5.1589012589541241E-006, 5.2724123375982313E-006, 5.3541463433315518E-006, 5.4319452973488174E-006, 5.5075421877648860E-006, 5.5823389327309226E-006, 5.6573264921061468E-006, 5.7331130970281323E-006, 5.8099970472481038E-006, 5.8880481629593514E-006, 5.9671797792140002E-006, 6.0472013522239294E-006, 6.1278581098764976E-006, 6.2088597847395375E-006, 6.2899004396107329E-006, 6.3706714335944688E-006, 6.4508692431604923E-006, 6.5301999362917668E-006, 6.6083812134502764E-006, 6.6851428103308001E-006, 6.7602267634381726E-006, 6.8333876646758279E-006, 6.9043918490528632E-006, 6.9730165925652158E-006, 7.0390495369602778E-006, 7.1022875296419982E-006, 7.1625361494518770E-006, 7.2196094152230139E-006, 7.2733296219871357E-006, 7.3235272517727152E-006, 7.3700411753602513E-006, 7.4127189379442331E-006, 7.4514162396539726E-006, 7.4859967502779660E-006, 7.5163315155465281E-006, 7.5422980716673646E-006, 7.5637799817002112E-006, 7.5806667531754469E-006, 7.5928538873846757E-006, 7.6002426334996788E-006, 7.6027391109787930E-006, 7.6002539777471323E-006, 7.5927031315267865E-006, 7.5800077173674855E-006, 7.5620942335531342E-006, 7.5388940064738222E-006, 7.5103432356458035E-006, 7.4763832308054822E-006, 7.4369606331940958E-006, 7.3920276158567506E-006, 7.3415420426939044E-006]
missing_phos_fortran = [-4.9957272742735606E-008, -5.0454978393331247E-008, -5.0941101536620295E-008, -5.1415617226168751E-008, -5.1878501590894929E-008, -5.2329731387397093E-008, -5.2769281594467596E-008, -5.3197126492745426E-008, -5.3613240345620858E-008, -5.4017597573092607E-008, -5.4410172662684859E-008, -5.4790942620494878E-008, -5.5159884497243434E-008, -5.5516977482923595E-008, -5.5862203192919519E-008, -5.6195544896863892E-008, -5.6516986058347810E-008, -5.6826509251134796E-008, -5.7124095800749839E-008, -5.7409726412429217E-008, -5.7683381895405517E-008, -5.7945044344576472E-008, -5.8194699347059203E-008, -5.8432336261192712E-008, -5.8657946830312958E-008, -5.8871523621955305E-008, -5.9073059987673925E-008, -5.9262550561483756E-008, -5.9439991637538070E-008, -5.9605381601734220E-008, -5.9758721264824076E-008, -5.9900014243587576E-008, -6.0029267257979809E-008, -6.0146490481855331E-008, -6.0251697757886490E-008, -6.0344906858732132E-008, -6.0426139615574857E-008, -6.0495422108097665E-008, -6.0552784769446553E-008, -6.0598262173500943E-008, -6.0631891993316091E-008, -6.0653714225856733E-008, -6.0663771011630989E-008, -6.0662107848117621E-008, -6.0648771317080291E-008, -6.0623808115403202E-008, -6.0587264780010741E-008, -6.0539187604952347E-008, -6.0479622617999197E-008, -6.0408615525790661E-008, -6.0326211915094299E-008, -6.0232457399497739E-008, -6.0127398009196575E-008, -6.0011080502201883E-008, -5.9883552898799576E-008, -5.9744864916028141E-008, -5.9595068604955875E-008, -5.9434218873656627E-008, -5.9262374070486163E-008, -5.9079595982522501E-008, -5.8885948398938528E-008, -5.8681495630005732E-008, -5.8466300301088997E-008, -5.8240423649811711E-008, -5.8003924619922639E-008, -5.7756858161680478E-008, -5.7499274105943326E-008, -5.7231217119693813E-008, -5.6952727518999675E-008, -5.6663842831077738E-008, -5.6364599441410035E-008, -5.6055033639654005E-008, -5.5735181068505236E-008, -5.5405075937610832E-008, -5.5064750044323107E-008, -5.4714231608778269E-008, -5.4353543796451900E-008, -5.3982703049885439E-008, -5.3601716968629524E-008, -5.3210581135030169E-008, -5.2809271428367839E-008, -5.2397740317232110E-008, -5.1975915758388121E-008, -5.1543698879219231E-008, -5.1100959464215354E-008, -5.0647529079524307E-008, -5.0183192548164031E-008, -4.9707679617433263E-008, -4.9220658397144796E-008, -4.8721731418016763E-008, -4.8210434433892037E-008, -4.7686244192853665E-008, -4.7148580909864540E-008, -4.6596826999513232E-008, -4.6030361751826386E-008, -4.5448613935129467E-008, -4.4851129231991955E-008, -4.4237637964938566E-008, -4.3608080016147110E-008, -4.2962498875351804E-008, -4.2214816050624246E-008, -4.1379919087224370E-008, -4.0552448942760530E-008, -3.9737029625135600E-008, -3.8937233709256498E-008, -3.8155645456921419E-008, -3.7393629390350867E-008, -3.6651566977113799E-008, -3.5929138791260594E-008, -3.5225573545795923E-008, -3.4539838715617690E-008, -3.3870778513554224E-008, -3.3217206846149662E-008, -3.2577966263296735E-008, -3.1951962423069169E-008, -3.1338182331029288E-008, -3.0735701320091658E-008, -3.0143682200266415E-008, -2.9561368906288448E-008, -2.8988078927561544E-008, -2.8423201503786064E-008, -2.7866187192945254E-008, -2.7316541616888558E-008, -2.6773821391643654E-008, -2.6237631238303883E-008, -2.5707620705229486E-008, -2.5183479677628861E-008, -2.4664934402888467E-008, -2.4151743650008856E-008, -2.3643695546599161E-008, -2.3140604482013675E-008, -2.2642308555846334E-008, -2.2148666941739484E-008, -2.1659557756607734E-008, -2.1174875938648599E-008, -2.0694531556254494E-008, -2.0218448013633852E-008, -1.9746560094171522E-008, -1.9278811395998950E-008, -1.8815152030744231E-008, -1.8355537792469954E-008, -1.7899931133984144E-008, -1.7448300828629705E-008, -1.7000622371276229E-008, -1.6556876819779322E-008, -1.6117050326588691E-008, -1.5681133986064885E-008, -1.5249123649609786E-008, -1.4821019731072175E-008, -1.4396826961287538E-008]
actual_uptake_po4_fortran = [2.3322032227480109E-006, 2.2159108475233858E-006, 2.0998711028945595E-006, 1.9842753398030213E-006, 1.8693113766360872E-006, 1.7551634056100194E-006, 1.6420118394842449E-006, 1.5300333072479151E-006, 1.4194006239194659E-006, 1.3102827304881900E-006, 1.2028446317255370E-006, 1.0972473799608120E-006, 9.9364776330898771E-007, 8.9219762015618414E-007, 7.9304387361634361E-007, 6.9632883884215754E-007, 6.0219048337636842E-007, 5.1076244145566155E-007, 4.2217397946064903E-007, 3.3654996770127697E-007, 2.5401084762582901E-007, 1.7467246886914725E-007, 9.8645668974234488E-008, 2.6036578892845950E-008, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 2.5010916228694582E-008, 1.0060151774625405E-007, 1.8012627835063843E-007, 2.6350673229612082E-007, 3.5065963302137560E-007, 4.4149695547807778E-007, 5.3592580033760004E-007, 6.3384839444999480E-007, 7.3516211394711313E-007, 8.3975950362100850E-007, 9.4752830529255471E-007, 1.0583515079658796E-006, 1.1721074286052046E-006, 1.2886698320185539E-006, 1.4079080881666861E-006, 1.5296872801970861E-006, 1.6538682266194152E-006, 1.7803070497872489E-006, 1.9088542237249201E-006, 2.0393541288224435E-006, 2.1716448100857106E-006, 2.3055581376139173E-006, 2.4409212041729019E-006, 2.5775565360750501E-006, 2.7152820982625041E-006, 2.8539115018167596E-006, 2.9932540151126812E-006, 3.1331135258763517E-006, 3.2732874696606639E-006, 3.4135668252029530E-006, 3.5537369395992472E-006, 3.6935801907979173E-006, 3.8328792119465435E-006, 3.9714204567534590E-006, 4.1089993985000854E-006, 4.2454265102210280E-006, 3.7727607321250627E-006, 3.3018245743879028E-006, 2.8579569014983826E-006, 2.4401796144253557E-006, 2.0475629391555822E-006, 1.6792179815269846E-006, 1.3342848106069366E-006, 1.0098245949112155E-006, 7.0703779471302943E-007, 4.2749351278915900E-007, 1.6987978583945385E-007, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000]

missing_assim_sum = np.zeros(len(assim_uptake_po4_fortran))
for i in range(0,len(assim_uptake_po4_fortran)):
    missing_assim_sum[i] = missing_phos_fortran[i] + assim_uptake_po4_fortran[i]