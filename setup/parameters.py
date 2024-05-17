# ----------------------------------------------------------------------------------------------------
# Superclass - PhytoParameters
# Subclasses - CoccolithophoreParameters, DiatomParameters

class PhytoParameters:
    """
    Description: Superclass to define phytoplankton parameters.
                 Additional parameters can be added for calcium carbonate (CaCO3) and
                 silica (Si) using the Coccolithophore and Diatom subclasses, respectively.
    """
    def __init__(self, constituents):
        
        self.chl_light_absorp = 0.
        self.chl_relax_rate = 0.
        self.chl_switch = 0
        self.max_chl_quota = 0.
        self.max_photo_rate = 0.

        self.spec_aff_nit = 0.
        self.nit_mult_factor = 0.
        self.max_nit_quota = 0.
        self.min_nit_quota = 0.
        self.opt_nit_quota = 0.
        self.half_sat_amm = 0.

        self.spec_aff_phos = 0.
        self.phos_mult_factor = 0.
        self.max_phos_quota = 0.
        self.min_phos_quota = 0.
        self.opt_phos_quota = 0.

        if "Fe" in constituents:
            self.spec_aff_iron = 0.
            self.iron_mult_factor = 0.
            self.min_iron_quota = 0.
            self.opt_iron_quota = 0.

        self.extra_lysis_rate = 0.
        self.half_sat_lysis = 0.
        self.half_sat_extra_lysis = 0.
        self.max_spec_lysis_rate = 0.

        self.max_light_util = 0.
        self.opt_par_irrad_ratio = 0.
        self.spec_resp_rate = 0.
        self.excr_prim_prod = 0.
        self.activity_resp_frac = 0.
        self.cuttof_temp_threshold = 0.
        self.nut_stress_threshold = 0.
        self.background_sink_rate = 0.
        self.bot_burial_vel = 0.
        self.max_sedi_rate = 0.

class CoccolithophoreParameters(PhytoParameters):
    """
    Description: Subclass for coccolithophore-specific parameters.
    """
    def __init__(self):
        self.caco3_frac = 0.
        self.caco3_prod = 0.
        self.caco3_quota = 0.

class DiatomParameters(PhytoParameters):
    """
    Description: Subclass for diatom-specific parameters.
    """
    def __init__(self):
        self.spec_aff_si = 0.
        self.half_sat_si = 0.
        self.half_sat_si_contois = 0.
        self.si_switch = 0
        self.min_si_quota = 0.
        self.opt_si_quota = 0.


# ----------------------------------------------------------------------------------------------------
# Superclass - ZooParameters
# Subclasses - MesoZooParameters, MicroZooParameters

class ZooParameters:
    """
    Definition: Superclass to define zooplankton parameters. 
                Additional meso- and micro-specific parameters provided in subclasses.
    """
    def __init__(self):
        self.assim_eff = 0. 
        self.excr_uptake = 0. 
        self.spec_growth_rate = 0. 
        self.spec_mort_rate = 0.
        self.spec_resp_rate = 0.
    
        self.max_nit_quota = 0.
        self.opt_nit_quota = 0.
        self.max_phos_quota = 0.
        self.opt_phos_quota = 0.
        self.half_sat_oxy = 0.

class MesoZooParameters(ZooParameters):
    """
    Definition: Subclass for mesozooplankton-specific parameters.
    """
    def __init__(self):
        self.nut_excr_rate = 0.
        self.dens_mort_rate = 0.
        self.dens_mort_rate_exp = 0.

class MicroZooParameters(ZooParameters):
    """
    Definition: Subclass for microzooplankton-specific parameters.
    """
    def __init__(self):
        self.feeding_threshold = 0.
        self.michaelis_constant = 0.
        self.mort_rate_oxy = 0.

