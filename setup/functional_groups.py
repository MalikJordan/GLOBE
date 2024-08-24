import numpy as np

# ----------------------------------------------------------------------------------------------------
# Living Organic Matter - Bacterioplankton, Phytoplankton, Zooplankton
class LivingOrganic:
    """
    Description: Superclass for the Living Organic functional groups.
                 Includes bacterioplankton, phytoplankton, and zooplankton.
                 Constituents included for each group are:
                 - Carbon (C)
                 - Nitrogen (N)
                 - Phosphorus (P)
    """
    def __init__(self, name, constituents, parameters, reactions, tracer):
        self.name = name
        self.conc_rate = {}
        self.parameters = {}
        self.sources = {}
        self.sinks = {}

        self.add_reactions(reactions, tracer)
        for const in constituents:
            self.add_constituent(const)
        for param in parameters:
            self.add_parameter(param,parameters[param])

    def add_constituent(self, constituent):
        self.conc_rate[constituent] = [0.,0.]

    def add_parameter(self, parameter, value):
        self.parameters[parameter] = value

    def add_reactions(self, reactions, tracer):
        for rate in reactions:
            if rate["consumed"] == tracer:
                self.sinks[rate["type"]] = rate["produced"]
            if rate["produced"] == tracer:
                self.sources[rate["type"]] = rate["consumed"]

# ----------------------------------------------------------------------------------------------------
# Superclass - Non-Living Organic Matter (Detritus) - Dissolved and Particulate Organic Matter
class NonLivingOrganic:
    """
    Description: Superclass for the Non-Living Organic functional group.
                 Includes dissolved and particulate organic matter.
    """
    def __init__(self, name, constituents, parameters, reactions, tracer):
        self.name = name
        self.conc_rate = {}
        self.parameters = {}
        self.sources = {}
        self.sinks = {}

        self.add_reactions(reactions, tracer)
        for const in constituents:
            self.add_constituent(const)
        for param in parameters:
            self.add_parameter(param,parameters[param])

    def add_constituent(self, constituent):
        self.conc_rate[constituent] = [0.,0.]

    def add_parameter(self, parameter, value):
        self.parameters[parameter] = value

    def add_reactions(self, reactions, tracer):
        for rate in reactions:
            if rate["consumed"] == tracer:
                self.sinks[rate["type"]] = rate["produced"]
            if rate["produced"] == tracer:
                self.sources[rate["type"]] = rate["consumed"]

# ----------------------------------------------------------------------------------------------------
# Class - Inorganic Tracers
class Inorganic:
    """
    Description: Class for inorganic tracers in the BGC model. Supported tracers include:
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
    def __init__(self, name, parameters, reactions, tracer):
        self.name = name
        self.constituents = [0.,0.]
        self.parameters = {}
        self.sources = {}
        self.sinks = {}

        self.add_reactions(reactions, tracer)
        for param in parameters:
            self.add_parameter(param, parameters[param])
        
    def add_parameter(self, parameter, value):
        self.parameters[parameter] = value
    
    def add_reactions(self, reactions, tracer):
        for rate in reactions:
            if rate["consumed"] == tracer:
                self.sinks[rate["type"]] = rate["produced"]
            if rate["produced"] == tracer:
                self.sources[rate["type"]] = rate["consumed"]