<!-- Description of phytoplankton species present in model. Phytoplankton species are distinguished by their 
equivalent spherical diameter (ESD) into three Phytoplankton Size Classes (PSCs). Example species provided
for each PSC.
-----------------------------
|    PSC    |   ESD (μm)    |
-----------------------------
|   Pico-   |     < 2       |   Cyanobacteria (Prochlorococcus, Synechococcus)
|   Nano-   |    2 - 20     |   Coccolithophores, Diazotrophs, NanoFlagellates
|   Micro-  |   20 - 200    |   Diatoms, DinoFlagellates
-----------------------------

The following contituents are included for phytoplanton groups.
Silica (Si) is specific to Diatoms and Calcium Carbonate (CaCO3) is specific to Coccolithophores.
-----------------------------------------
|     Constituent   |      Units        |
-----------------------------------------
|       Carbon      |   mmol C / m3     |
|      Nitrogen     |   mmol N / m3     |
|     Phosphorus    |   mmol P / m3     |
|   Chlorophyll-a   |   mg Chl / m3     |
|        Iron       |   μmol Fe / m3    |
-----------------------------------------
|       Silica      |   mmol Si / m3    |   Diatoms Only
| Calcium Carbonate |   mmol C / m3     |   Coccolithophores Only
-----------------------------------------

The following rates/reactions are considered for phytoplankton groups.
-----------------------------------------------------------------------------------------
|        Rate/Reaction          |   Source/Sink     |      Constituents Affected        |
-----------------------------------------------------------------------------------------
|        Calcification          |     Source        |                   , CaCO3         |
|   Gross Primary Production    |     Source        |   C                               |
|          Exudation            |      Sink         |   C                               |
|            Lysis              |      Sink         |   C, N, P,      Fe, CaCO3, Si     |
|       Nutrient Uptake         |     Source        |      N, P,      Fe,        Si     |
|          Predation            |      Sink         |   C, N, P, Chl, Fe, CaCO3, Si     |
|         Respiration           |      Sink         |   C                               |
|          Synthesis            |     Source        |            Chl                    |
----------------------------------------------------------------------------------------- 

The following parameters are used in rate equations for phytoplankton groups. Calcium Carbonate (CaCO3)
parameters are only used for Coccolithophores. Silica (Si) parameters are only used for Diatoms.
-----------------------------------------------------------------------------------------------------------------------------------------------------
|   Name                    |   Description                                                                                                         |      
----------------------------------------------------------------------------------------------------------------------------------------------------- 
|   chl_light_absorp        |   Chl-specific light absorption coefficient                                                                           |
|   chl_relax_rate          |   Maximum specific photosynthetic rate                                                                                |
|   chl_switch              |   Switch for the type of dissolved organic carbon ecretion (Consistent with bacteria parameterization [bac_version])  |
|   max_chl_quota           |   Maximum Chl:C quota                                                                                                 |
|   max_photo_rate          |   Maximum specific photosynthetic rate                                                                                |
-----------------------------------------------------------------------------------------------------------------------------------------------------
|   spec_aff_nit            |   Specific affinity constant for nitrogen                                                                             |
|   nit_mult_factor         |   Multiplication factor for luxury storage of nitrogen                                                                |
|   max_nit_quota           |   Maximum nitrogen quota                                                                                              |
|   min_nit_quota           |   Minimum nitrogen quota                                                                                              |
|   opt_nit_quota           |   Optimal/Reference nitrogen quota                                                                                    |
|   half_sat_amm            |   Half saturation constant for ammonium uptake preference over ammonium                                               |
-----------------------------------------------------------------------------------------------------------------------------------------------------
|   spec_aff_phos           |   Specific affinity constant for phosphorus                                                                           |
|   phos_mult_factor        |   Multiplication factor for luxury storage of phosophorus                                                             |
|   max_phos_quota          |   Maximum phosphorus quota                                                                                            |
|   min_phos_quota          |   Minimum phosphorus quota                                                                                            |
|   opt_phos_quota          |   Optimal/Reference phosphorus quota                                                                                  |
-----------------------------------------------------------------------------------------------------------------------------------------------------
|   spec_aff_iron           |   Specific affinity constant for iron                                                                                 |
|   iron_mult_factor        |   Multiplication factor for luxury storage of iron                                                                    |
|   min_iron_quota          |   Minimum quotum for Fe:C                                                                                             |
|   opt_iron_quota          |   Optimal/Reference Fe:C ratio                                                                                        |
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
|   caco3_frac              |   Fraction of phytoplankton organic matter associated with CaCO3                                                      |   Coccolithophores Only
|   caco3_prod              |   Calcification rate                                                                                                  |
|   caco3_quota             |   CaCO3:C ratio                                                                                                       |
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
|   spec_aff_si             |   Specific affinity constant for silica                                                                               |   Diatoms Only
|   half_sat_si             |   Half saturation constant for Si-limitation                                                                          |
|   half_sat_si_contois     |   Variable half saturation constant for Si-limitation of Contois                                                      |
|   si_switch               |   Switch for parameterixation of silicate limitation (1=external, 2=internal)                                         |
|   min_si_quota            |   Minimum SI:C ratio in silicifiers                                                                                   |
|   opt_si_quota            |   Optimal/Reference Si:C ratio in slicifiers                                                                          |
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
|   extra_lysis_rate        |   Extra lysis rate                                                                                                    |
|   half_sat_lysis          |   Half saturation constant for nutrient stress lysis (Nutrient stress threshold)                                      |
|   half_sat_extra_lysis    |   Half saturation constant for extra lysis                                                                            |
|   max_spec_lysis_rate     |   Maximum specific nutrient-stress lysis rate                                                                         |
-----------------------------------------------------------------------------------------------------------------------------------------------------
|   max_light_util          |   Maximum light utilization coefficient                                                                               |
|   opt_par_irrad_ratio     |   Optimal value of E_PAR / E_K                                                                                        |
|   spec_resp_rate          |   Basal specific respiration rate                                                                                     |
|   excr_prim_prod          |   Excreted fraction of primary production                                                                             |
|   activity_resp_frac      |   Activity respiration fraction                                                                                       |
|   cuttof_temp_threshold   |   Cut-off threshold for temperature regulating factor                                                                 |
|   nut_stress_threshold    |   Nutrient stress threshold for sinking                                                                               |
|   background_sink_rate    |   Additional background sinking rate                                                                                  |
|   bot_burial_vel          |   Bottom burial velocity for plankton                                                                                 |
|   max_sedi_rate           |   Maximum sedimentation rate                                                                                          |
-----------------------------------------------------------------------------------------------------------------------------------------------------
-->
