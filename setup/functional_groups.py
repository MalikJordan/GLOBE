import numpy as np

# ----------------------------------------------------------------------------------------------------
# Superclass - LivingOrganic
# Subclasses - Bacterioplankton (Non-Phototrophic), Phytoplankton, Zooplankton

class LivingOrganic:
    """
    Description: Superclass for the Living Organic functional groups.
                 Includes bacterioplankton, phytoplankton, and zooplankton.
                 Constituents included for each group are:
                 - Carbon (C)
                 - Nitrogen (N)
                 - Phosphorus (P)
    """
    def __init__(self, name, num_boxes):
        self.name = name

        if num_boxes > 1:
            self.C = np.zeros(num_boxes,dtype=float)
            self.dC_dt = np.zeros(num_boxes,dtype=float)

            self.N = np.zeros(num_boxes,dtype=float)
            self.dN_dt = np.zeros(num_boxes,dtype=float)

            self.P = np.zeros(num_boxes,dtype=float)
            self.dP_dt = np.zeros(num_boxes,dtype=float)
        else:
            self.C = 0.
            self.dC_dt = 0.

            self.N = 0.
            self.dN_dt = 0.

            self.P = 0.
            self.dP_dt = 0.

class Bacterioplankton(LivingOrganic):
    """
    Description: Creates subclass for bacterioplankton groups in the BGC model.
                 No additional constituents are added to those already included in LivingOrganic Superclass.
    """
    def __init__(self, name, num_boxes):
        super().__init__(name, num_boxes)

class Phytoplankton(LivingOrganic):
    """
    Description: Creates subclass for phytoplankton groups in the BGC model.
                 Chlorophyll-a is an added constituent for all phytoplankton groups.
                 Additional constituents can be added for:
                 - Iron (Fe)
                 - Silica (Si)
    """
    def __init__(self, name, constituents, num_boxes):
        super().__init__(name, num_boxes)

        if num_boxes > 1:
            self.Chl = np.zeros(num_boxes,dtype=float)
            self.dChl_dt = np.zeros(num_boxes,dtype=float)

            if "Fe" in constituents:
                self.Fe = np.zeros(num_boxes,dtype=float)
                self.dFe_dt = np.zeros(num_boxes,dtype=float)
            if "Si" in constituents:
                self.Si = np.zeros(num_boxes,dtype=float)
                self.dSi_dt = np.zeros(num_boxes,dtype=float)
        else:
            self.Chl = 0.
            self.dChl_dt = 0.

            if "Fe" in constituents:
                self.Fe = 0.
                self.dFe_dt = 0.
            if "Si" in constituents:
                self.Si = 0.
                self.dSi_dt = 0.

class Zooplankton(LivingOrganic):
    """
    Description: Creates subclass for zooplankton groups in the BGC model.
                 No additional constituents are added to those already included in LivingOrganic Superclass.
    """
    def __init__(self, name, num_boxes):
        super().__init__(name, num_boxes)

# ----------------------------------------------------------------------------------------------------
# Superclass - NonLivingOrganic
# Subclasses - OrganicMatter (includes dissolved and particulate)

class NonLivingOrganic:
    """
    Description: Superclass for the Non-Living Organic functional group.
                 Includes dissolved and particulate organic matter.
    """
    def __init__(self, name):
        self.name = name

class OrganicMatter(NonLivingOrganic):
    """
    Description: Creates subclasses for dissolved and particulate inorganic matter.
                 Constituents which may be included in nutrient and detrital pools include:
                 - Carbon (C)
                 - Iron (Fe)
                 - Nitrogen (N)
                 - Phosphorus (P)
                 - Silicate (Si)
    """
    def __init__(self, name, constituents, num_boxes):
        super().__init__(name)

        if num_boxes > 1:
            if "C" in constituents:
                self.C = np.zeros(num_boxes,dtype=float)
                self.dC_dt = np.zeros(num_boxes,dtype=float)
            if "N" in constituents:
                self.N = np.zeros(num_boxes,dtype=float)
                self.dN_dt = np.zeros(num_boxes,dtype=float)
            if "P" in constituents:
                self.P = np.zeros(num_boxes,dtype=float)
                self.dP_dt = np.zeros(num_boxes,dtype=float)
            if "Si" in constituents:
                self.Si = np.zeros(num_boxes,dtype=float)
                self.dSi_dt = np.zeros(num_boxes,dtype=float)
            if "Fe" in constituents:
                self.Fe = np.zeros(num_boxes,dtype=float)
                self.dFe_dt = np.zeros(num_boxes,dtype=float)
        else:
            if "C" in constituents:
                self.C = 0.
                self.dC_dt = 0.
            if "N" in constituents:
                self.N = 0.
                self.dN_dt = 0.
            if "P" in constituents:
                self.P = 0.
                self.dP_dt = 0.
            if "Si" in constituents:
                self.Si = 0.
                self.dSi_dt = 0.
            if "Fe" in constituents:
                self.Fe = 0.
                self.dFe_dt = 0.

# ----------------------------------------------------------------------------------------------------
# Class - Inorganic
class Inorganic:
    """
    Description: Creates a class of inorganic tracers in the BGC model. Supported tracers include:
                 - Ammonium (NH4)
                 - Calcium Carbonate (CaCO3)
                 - Dissolved Inorganic Carbon (DIC)
                 - Iron (Fe)
                 - Nitrate (NO3)
                 - Nitrogen Sink (N2)
                 - Oxygen (O2)
                 - Phosphate (PO4)
                 - Silicate (SiO2)
                 - Sulfide (S)
                 - Total Alkalinity (TALK)
    """
    def __init__(self, tracers, num_boxes):
        if num_boxes > 1:
            if "CaCO3" in tracers:
                self.CaCO3 = np.zeros(num_boxes,dtype=float)
                self.dCaCO3_dt = np.zeros(num_boxes,dtype=float)
            if "DIC" in tracers:
                self.DIC = np.zeros(num_boxes,dtype=float)
                self.dDIC_dt = np.zeros(num_boxes,dtype=float)
            if "Fe" in tracers:
                self.Fe = np.zeros(num_boxes,dtype=float)
                self.dFe_dt = np.zeros(num_boxes,dtype=float)
            if "N2" in tracers:
                self.N2 = np.zeros(num_boxes,dtype=float)
                self.dN2_dt = np.zeros(num_boxes,dtype=float)
            if "NH4" in tracers:
                self.NH4 = np.zeros(num_boxes,dtype=float)
                self.dNH4_dt = np.zeros(num_boxes,dtype=float)
            if "NO3" in tracers:
                self.NO3 = np.zeros(num_boxes,dtype=float)
                self.dNO3_dt = np.zeros(num_boxes,dtype=float)
            if "O2" in tracers:
                self.O2 = np.zeros(num_boxes,dtype=float)
                self.dO2_dt = np.zeros(num_boxes,dtype=float)
            if "PO4" in tracers:
                self.PO4 = np.zeros(num_boxes,dtype=float)
                self.dPO4_dt = np.zeros(num_boxes,dtype=float)
            if "S" in tracers:
                self.S = np.zeros(num_boxes,dtype=float)
                self.dS_dt = np.zeros(num_boxes,dtype=float)
            if "SiSO2" in tracers:
                self.SiO2 = np.zeros(num_boxes,dtype=float)
                self.dSiO2_dt = np.zeros(num_boxes,dtype=float)
            if "TALK" in tracers:
                self.TALK = np.zeros(num_boxes,dtype=float)
                self.dTALK_dt = np.zeros(num_boxes,dtype=float)
        else:
            if "CaCO3" in tracers:
                self.CaCO3 = 0.
                self.dCaCO3_dt = 0.
            if "DIC" in tracers:
                self.DIC = 0.
                self.dDIC_dt = 0.
            if "Fe" in tracers:
                self.Fe = 0.
                self.dFe_dt = 0.
            if "N2" in tracers:
                self.N2 = 0.
                self.dN2_dt = 0.
            if "NH4" in tracers:
                self.NH4 = 0.
                self.dNH4_dt = 0.
            if "NO3" in tracers:
                self.NO3 = 0.
                self.dNO3_dt = 0.
            if "O2" in tracers:
                self.O2 = 0.
                self.dO2_dt = 0.
            if "PO4" in tracers:
                self.PO4 = 0.
                self.dPO4_dt = 0.
            if "S" in tracers:
                self.S = 0.
                self.dS_dt = 0.
            if "SiSO2" in tracers:
                self.SiO2 = 0.
                self.dSiO2_dt = 0.
            if "TALK" in tracers:
                self.TALK = 0.
                self.dTALK_dt = 0.