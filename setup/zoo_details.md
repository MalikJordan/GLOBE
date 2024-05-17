<!-- Description of zooplankton species present in model. Zooplankton species are distinguished by their 
equivalent spherical diameter (ESD) into two Zooplankton Size Classes (ZSCs). Example species provided
for each ZSC.
---------------------------------
|    ZSC    |        ESD        |
---------------------------------
|   Micro-  |   20 - 200 μm     |   Ciliates, DinoFlagellates
|   Meso-   |   0.2 - 20 mm     |   Copepods, Decapod Larvae, Gastropod Larvae
---------------------------------

The following rates/reactions are considered for zooplankton groups. 
Mesozooplankton groups can prey on phytoplankton, microzooplankton, and mesozooplankton.
Microzooplankton groups can prey on bacterioplankton, phytoplankton, and microzooplankton.
-----------------------------------------------------------------------------
|   Rate/Reaction   |   Source/Sink     |      Constituents Affected        |
-----------------------------------------------------------------------------
|     Predation     |   Source/Sink     |              C, N, P              |
|      Release      |      Sink         |              C, N, P              |
|    Respiration    |      Sink         |              C                    |
-----------------------------------------------------------------------------

The following parameters are used in rate equations for zooplankton groups.
-----------------------------------------------------------------------------------------
|   Name                    |   Description                                             |
-----------------------------------------------------------------------------------------
|   assim_eff               |   Assimilation efficiency                                 |
|   excr_uptake             |   Excreted fraction of uptake                             |
|   spec_growth_rate        |   Potential specific growth rate                          |
|   spec_mort_rate          |   Specific mortality rate                                 |
|   spec_resp_rate          |   Basal respiration rate                                  |
|   max_nit_quota           |   Maximum nitrogen quota                                  |
|   opt_nit_quota           |   Optimal nitrogen quota                                  |
|   max_phos_quota          |   Maximum phosphorus quota                                |
|   opt_phos_quota          |   Optimal phosphorus quota                                |
|   half_sat_oxy            |   Half saturation value for oxygen                        |
-----------------------------------------------------------------------------------------------------------------------
|   feeding_threshold       |   Feeding threshold                                       |   Microzooplankton Only
|   michaelis_constant      |   Michaelis constant for total food ingestion             |
|   mort_rate_oxy           |   Mortality rate due to oxygen limitation                 |
-----------------------------------------------------------------------------------------------------------------------
|   dens_mort_rate          |   Density-dependent specific mortality rate               |   Mesozooplankton Only
|   dens_mort_rate_exp      |   Exponent for density-dependent mortality rate           |
|   nut_excr_rate           |   Specific rate of nutrients and carbon excretion         |
-----------------------------------------------------------------------------------------
-->