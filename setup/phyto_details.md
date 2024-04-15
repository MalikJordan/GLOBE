<!-- Initialization of phytoplankton species present in model. Phytoplankton species are distinguished by their 
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
----------------------------------------------------------------------------------------- -->
