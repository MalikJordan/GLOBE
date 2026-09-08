import numpy as np
from numba import njit, types
from numba.types import float64, unicode_type
from numba.typed import Dict, List
from pom.calculations import temperature_and_salinity_profiles
from functions.bgc_rate_eqns import bgc_rate_eqns
import os
from pom.compare_values import compare_values
from tests.bfm17 import check_conc, check_rates
# from tests.bfm17 import check_after_vdiff_SOS
from tests.bfm56 import check_after_vdiff_SOS
np.set_printoptions(precision=20)

phytoc0 = [0.13600693026060057,
0.13788003350628575,
0.13783314465202093,
0.14077789534980695,
0.1406810805896549,
0.14682565519044186,
0.1469538277508386,
0.15355533110108344,
0.15380860337581495,
0.15782714350506166,
0.15823479872248547,
0.1603634138729713,
0.16088357936184067,
0.16229641730321456,
0.16290433707771607,
0.1640265756317987,
0.1646512411912123,
0.1656892933940167,
0.16631613267471046,
0.1672801529436927,
0.16791557987374542,
0.16881654457000914,
0.16945076008004398,
0.17030422937999282,
0.17093152355074928,
0.1717461382009736,
0.17236375721963207,
0.1731440730732346,
0.1737537447907273,
0.17450355689311903,
0.17510471623443227,
0.1758270206516259,
0.17642106478320435,
0.17711808642920107,
0.177705748745344,
0.1783798852086758,
0.17896198303280988,
0.17961509006864942,
0.1801919679459737,
0.18082550828448518,
0.18139692518735268,
0.182012595189088,
0.18257749847798416,
0.18317799559609163,
0.18373411837356493,
0.18432338606561416,
0.18486732106459608,
0.1854495871781066,
0.1859787882754807,
0.1865556972379106,
0.18707352725578905,
0.1876417647664717,
0.1881523512368167,
0.1887092230186375,
0.18921355314393296,
0.1897590418370027,
0.19025716311759905,
0.190791741597451,
0.19128346616931985,
0.19180778594882614,
0.19229288917344098,
0.19280768371838275,
0.19328581852583587,
0.19379192274072263,
0.194262651965612,
0.1947609137088863,
0.19522384904008236,
0.1957149463426617,
0.19616996978725082,
0.19665421630946184,
0.1971016047149702,
0.19757892343012395,
0.19801920670128392,
0.19848934083368955,
0.19892303144177167,
0.19938579787071883,
0.1998132372438923,
0.2002686325704546,
0.20068998025642543,
0.20113816460355233,
0.20155344655522792,
0.2019946876590887,
0.2024038543466883,
0.20283846568374028,
0.20324145270409028,
0.20366973081010561,
0.2040665180679751,
0.20448868860114627,
0.20487933962380078,
0.20529553244377985,
0.20568019346555458,
0.20609045976182916,
0.20646932114801636,
0.2068736799088542,
0.20724692789170784,
0.20764541145421228,
0.20801319731004164,
0.20840587461636734,
0.20876830855703515,
0.2091552847541635]

dphytoc0_dt = [2.4755744070392944e-06 ,
2.2818396418391657e-06 ,
2.2830341677615567e-06 ,
2.2686918872596887e-06 ,
2.2631694970350074e-06 ,
2.3088564741477703e-06 ,
2.3005284577016885e-06 ,
2.363951069053597e-06 ,
2.353414212408722e-06 ,
2.3787004279543424e-06 ,
2.3696736578622503e-06 ,
2.3648220235825352e-06 ,
2.356605711214298e-06 ,
2.341537671875432e-06 ,
2.334906171846532e-06 ,
2.3156927041796382e-06 ,
2.3073953662189327e-06 ,
2.290994489265538e-06 ,
2.281223684531372e-06 ,
2.264658357630849e-06 ,
2.2558781161608068e-06 ,
2.2389947710257535e-06 ,
2.2296407128710994e-06 ,
2.212967684752981e-06 ,
2.2030348032225503e-06 ,
2.1868218986392617e-06 ,
2.17655538456412e-06 ,
2.160814793884954e-06 ,
2.1501336671641675e-06 ,
2.1348928277316772e-06 ,
2.1239528780967814e-06 ,
2.1091385052853723e-06 ,
2.0979410461269117e-06 ,
2.0835443499021336e-06 ,
2.0721490724615306e-06 ,
2.0581576345582797e-06 ,
2.046618290433258e-06 ,
2.0330190124558263e-06 ,
2.02137422716915e-06 ,
2.008139435767859e-06 ,
1.9964271632435622e-06 ,
1.9835094894301315e-06 ,
1.9717970978552274e-06 ,
1.9591176290273324e-06 ,
1.9475207547016606e-06 ,
1.934966935424715e-06 ,
1.923634666677869e-06 ,
1.9110932665015686e-06 ,
1.9001260139055388e-06 ,
1.8875769911782143e-06 ,
1.8768559057773658e-06 ,
1.8644408980059532e-06 ,
1.8538214454643014e-06 ,
1.8416445345551042e-06 ,
1.8310944160824061e-06 ,
1.8191601048320136e-06 ,
1.8086824106099311e-06 ,
1.7969771336869415e-06 ,
1.7865851763037274e-06 ,
1.775090849882436e-06 ,
1.7648014729909736e-06 ,
1.7534994389117582e-06 ,
1.7433346005846775e-06 ,
1.7322045376119722e-06 ,
1.7221883409255328e-06 ,
1.7112102560978293e-06 ,
1.701362676678202e-06 ,
1.6905217165558564e-06 ,
1.6808506665341032e-06 ,
1.6701423414158121e-06 ,
1.6606402108366028e-06 ,
1.6500706040408564e-06 ,
1.6407201770962077e-06 ,
1.6302992212878232e-06 ,
1.6210836202099393e-06 ,
1.610817583308048e-06 ,
1.6017259202439877e-06 ,
1.5916148079521084e-06 ,
1.5826428104608803e-06 ,
1.5726816051006693e-06 ,
1.5638300785355485e-06 ,
1.554011039992666e-06 ,
1.545283838110929e-06 ,
1.5355987342470363e-06 ,
1.5270005711428142e-06 ,
1.5174426501061089e-06 ,
1.508977003867226e-06 ,
1.4995424169148317e-06 ,
1.4912101731463094e-06 ,
1.4818982941662773e-06 ,
1.4736977139891325e-06 ,
1.4645101050009422e-06 ,
1.4564379194967568e-06 ,
1.4473765001009257e-06 ,
1.4394292115735164e-06 ,
1.4304946653826841e-06 ,
1.4226692864403713e-06 ,
1.4138603737862277e-06 ,
1.406154494834026e-06 ,
1.3974682583182601e-06 ]

bfm_rates = [[
-0.0005622903659822886 ,
3.138844740080756e-08 ,
-2.5068995748224415e-09 ,
4.971958393648962e-07 ,
-5.530190167704545e-07 ,
2.4755744070392944e-06 ,
8.693495065911221e-09 ,
2.1018217412351313e-11 ,
7.841456670377728e-09 ,
-2.469816571466231e-06 ,
-3.112039324134038e-08 ,
-1.941813466042645e-09 ,
-1.077555360146571e-05 ,
-1.3747579582383005e-07 ,
-8.572681860471073e-09 ,
-2.6563039122260403e-05 ,
-3.3478624579081466e-07 ,
-2.0894970291706202e-08 ,
],
[
-0.00035109739807150277 ,
3.109085102529036e-08 ,
-2.8540896720075545e-09 ,
4.926844208668044e-07 ,
-5.526980552011531e-07 ,
2.2818396418391657e-06 ,
8.785869338179739e-09 ,
2.12349479542908e-11 ,
6.474147947478278e-09 ,
-2.3193503374292606e-06 ,
-2.922449278283736e-08 ,
-1.8235153063635402e-09 ,
-1.0705657505945826e-05 ,
-1.3739815939445633e-07 ,
-8.567828414210396e-09 ,
-2.634133481663457e-05 ,
-3.3199354835568285e-07 ,
-2.0720742252670714e-08 ,
],
[
-0.000357518512831956 ,
3.1083727921711784e-08 ,
-2.8487523083910906e-09 ,
4.92569988105756e-07 ,
-5.527043849691009e-07 ,
2.2830341677615567e-06 ,
8.780637511080873e-09 ,
2.122552455604397e-11 ,
6.484825734287621e-09 ,
-2.3193062047183394e-06 ,
-2.9223936255214085e-08 ,
-1.8234805809859055e-09 ,
-1.070738189314396e-05 ,
-1.3739226663777754e-07 ,
-8.567467392420082e-09 ,
-2.633278399106588e-05 ,
-3.3188567041545415e-07 ,
-2.0714005472861835e-08 ,
],
[
-0.000358612087963394 ,
3.120442874088037e-08 ,
-2.7899852340117448e-09 ,
4.947839526754338e-07 ,
-5.525108502915307e-07 ,
2.2686918872596887e-06 ,
8.423694217034529e-09 ,
2.007904878582083e-11 ,
6.27850515883561e-09 ,
-2.3412800793998713e-06 ,
-2.9500842908494207e-08 ,
-1.8407586629757163e-09 ,
-1.0617790596481304e-05 ,
-1.3713058750561342e-07 ,
-8.55109877297591e-09 ,
-2.6483520923865232e-05 ,
-3.33786231244349e-07 ,
-2.083265035371457e-08 ,
],
[
-0.00035807789610300803 ,
3.114582589122636e-08 ,
-2.8049201977847217e-09 ,
4.937871377266315e-07 ,
-5.52469082986163e-07 ,
2.2631694970350074e-06 ,
8.49986147904011e-09 ,
2.0306996086612133e-11 ,
6.255414927657017e-09 ,
-2.3379029858387204e-06 ,
-2.9458289500257597e-08 ,
-1.8381034661908385e-09 ,
-1.0615598833472438e-05 ,
-1.3714987465242234e-07 ,
-8.552311076545056e-09 ,
-2.6411124454273308e-05 ,
-3.3287391485520697e-07 ,
-2.077571834457708e-08 ,
],
[
-0.0003752665818835991 ,
3.1639713468201996e-08 ,
-2.383108130799395e-09 ,
5.02857139213085e-07 ,
-5.523669498893376e-07 ,
2.3088564741477703e-06 ,
6.844237006471142e-09 ,
1.5376006624380294e-11 ,
6.230174264538123e-09 ,
-2.4166912927924013e-06 ,
-3.045111170292951e-08 ,
-1.9000523845737113e-09 ,
-1.0470164510265753e-05 ,
-1.361531124900626e-07 ,
-8.489953065341617e-09 ,
-2.703309113931567e-05 ,
-3.407140438957646e-07 ,
-2.1265084024911053e-08 ,
],
[
-0.0003769640383178212 ,
3.164107950483612e-08 ,
-2.394453321353933e-09 ,
5.028369572434554e-07 ,
-5.523275669843113e-07 ,
2.3005284577016885e-06 ,
6.900295703560546e-09 ,
1.554153531681241e-11 ,
6.183281681097236e-09 ,
-2.4142880538089856e-06 ,
-3.04208314574219e-08 ,
-1.898162993159106e-09 ,
-1.0463120994259812e-05 ,
-1.361991325931837e-07 ,
-8.49283284923432e-09 ,
-2.7033791869320973e-05 ,
-3.407228355750564e-07 ,
-2.1265625197759508e-08 ,
],
[
-0.00039426791592601115 ,
3.213689957938693e-08 ,
-1.7606956921209657e-09 ,
5.120687661724612e-07 ,
-5.522340643856784e-07 ,
2.363951069053597e-06 ,
4.8837900132310866e-09 ,
9.471377168283964e-12 ,
6.226572716407266e-09 ,
-2.5051757339126082e-06 ,
-3.1566122718635774e-08 ,
-1.9696254900908206e-09 ,
-1.029489827543615e-05 ,
-1.347943431002449e-07 ,
-8.404983452434143e-09 ,
-2.7677025659191022e-05 ,
-3.488313946746907e-07 ,
-2.177176201403025e-08 ,
],
[
-0.00039591372259588305 ,
3.214705416004671e-08 ,
-1.763100246848316e-09 ,
5.121958669953381e-07 ,
-5.521836338409e-07 ,
2.353414212408722e-06 ,
4.9240671663099915e-09 ,
9.61111827018977e-12 ,
6.162951586769304e-09 ,
-2.503034405463463e-06 ,
-3.153914414752849e-08 ,
-1.967942112046762e-09 ,
-1.0283461427933966e-05 ,
-1.3484226059473362e-07 ,
-8.407981292188216e-09 ,
-2.7688454572758568e-05 ,
-3.489754291725377e-07 ,
-2.178074187408192e-08 ,
],
[
-0.00040481056943918913 ,
3.2392535722272814e-08 ,
-1.3660691006337756e-09 ,
5.167762870228568e-07 ,
-5.520622557958231e-07 ,
2.3787004279543424e-06 ,
3.826258649588383e-09 ,
6.203050145371472e-12 ,
6.115953300262629e-09 ,
-2.5525230300539776e-06 ,
-3.216276784759683e-08 ,
-2.0068542247854025e-09 ,
-1.0165753807648383e-05 ,
-1.339664911460235e-07 ,
-8.353230248409845e-09 ,
-2.801619563359032e-05 ,
-3.53107217578191e-07 ,
-2.2038654299222934e-08 ,
]]


globe_rates = np.zeros((100,18))

# def pom_bgc_1d(iter, base_element, light_attenuation_water, temperature, salinity, density, inorganic_suspended_matter, shortwave_radiation,
#                 w_eddy_velocity, w_gen, wind_speed_zonal, wind_speed_meridional, dif_trac,
#                 dt2, num_layers, vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, column_depth, 
#                 nrt_o2, nrt_po4, nrt_no3, nrt_nh4, o2b, no3b, ponb_grad, po4b,
#                 smoth, umolbgc, nbcbgc, ntp, water_specific_heat_times_density, 
#                 concentration, sinking, tracer_map, tracer_type, tracers):
    
#     # Extract current concentration and initialize rate of change array
#     conc = concentration[...,iter]
#     d_dt = np.zeros_like(conc,dtype=np.float64)

def pom_bgc_1d(iter, configuration, base_element, light_attenuation_water, temperature, salinity, density, inorganic_suspended_matter, shortwave_radiation,
                w_eddy_velocity, w_gen, wind_speed_zonal, wind_speed_meridional, dif_trac,
                dt2, num_layers, vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, column_depth, 
                nrt_o2, nrt_po4, nrt_no3, nrt_nh4, o2b, no3b, ponb_grad, po4b,
                smoth, umolbgc, nbcbgc, ntp, water_specific_heat_times_density, 
                conc_bwd, conc_cur, sinking, tracer_map, tracer_type, tracers):
    
    # Initialize rate of change array
    d_dt = np.zeros_like(conc_cur,dtype=np.float64)

    # # Reset sinking matrix to initial value
    # if iter == 0:   reset_sinking = sinking.copy()
    # else:   sinking = reset_sinking.copy()

    # if iter < 5:
    #     conc_name = f"conc_iter{iter:01d}.npy"
    #     load_conc = np.load(os.getcwd() + "/tests/bfm56/check_conc/" + conc_name, allow_pickle=True)
    #     # delta_conc = conc - load_conc
    #     delta_conc = conc_cur - load_conc

    #     x = 1

    # if iter > 24 and iter < 30:
    #     # concentrations are fine at iter==27, incorrect starting at iter==28
    #     # first group with wrong concentrations is phytoplankton (all constituents)
    #     # happens before rate calculations at iter==28
    #     # rates before this step are fine (match at iter==27 inside of vdiff_SOS)
    #     # check phyto sedimentation
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter24to29[iter-25]
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter24to29[iter-25]
    #     x = 1

    # Physical variables for bgc rate equations
    temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind = bgc_physical(temperature, salinity, density, inorganic_suspended_matter, vertical_spacing, column_depth, shortwave_radiation, water_specific_heat_times_density, w_eddy_velocity, w_gen, wind_speed_zonal, wind_speed_meridional)

    # Calculate rate of change
    # d_dt = bgc_rate_eqns(iter, base_element, conc, d_dt, light_attenuation_water, temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind, tracer_map, tracer_type, tracers, sinking)
    # d_dt = bgc_rate_eqns(iter, base_element, conc_cur, d_dt, light_attenuation_water, temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind, tracer_map, tracer_type, tracers, sinking)
    d_dt = bgc_rate_eqns(iter, configuration, base_element, conc_cur, d_dt, light_attenuation_water, temp, sal, dens, z, dz, surface_PAR, wind, tracer_map, tracer_type, tracers, sinking)
    
    if "o2" in tracers: d_o2surf = tracers["o2"].surf_flux
    else:   d_o2surf = np.float64(0.)

    # if iter < 5:
    #     rates_name = f"rates_iter{iter:01d}.npy"
    #     # phys_name = f"phys_iter{iter:01d}.npz"

    #     load_rates = np.load(os.getcwd() + "/tests/bfm56/check_rates/" + rates_name, allow_pickle=True)
    #     # load_phys = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + phys_name, allow_pickle=True)

    #     delta_rates = d_dt - load_rates
    # #     delta_temp = temp - load_phys["temp"]
    # #     delta_sal = sal - load_phys["sal"]
    # #     delta_dens = dens - load_phys["dens"]
    # #     delta_ism = ism - load_phys["ism"]
    # #     delta_z = z - load_phys["z"]
    # #     delta_dz = dz - load_phys["dz"]
    # #     delta_surface_PAR = surface_PAR - load_phys["surface_PAR"]
    # #     delta_weddy = weddy - load_phys["weddy"]
    # #     delta_wgen = wgen - load_phys["wgen"]
    # #     delta_wind = wind - load_phys["wind"]

    #     x = 1

    # if iter == 119:
    #     conc_name = f"conc_iter{iter:03d}.npy"
    #     rates_name = f"rates_iter{iter:03d}.npy"
    #     phys_name = f"phys_iter{iter:03d}.npz"

    #     load_conc = np.load(os.getcwd() + "/tests/bfm17/check_conc/" + conc_name, allow_pickle=True)
    #     load_rates = np.load(os.getcwd() + "/tests/bfm17/check_rates/" + rates_name, allow_pickle=True)
    #     load_phys = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + phys_name, allow_pickle=True)

    #     delta_conc = conc - load_conc
    #     delta_rates = d_dt - load_rates
    #     delta_temp = temp - load_phys["temp"]
    #     delta_sal = sal - load_phys["sal"]
    #     delta_dens = dens - load_phys["dens"]
    #     delta_ism = ism - load_phys["ism"]
    #     delta_z = z - load_phys["z"]
    #     delta_dz = dz - load_phys["dz"]
    #     delta_surface_PAR = surface_PAR - load_phys["surface_PAR"]
    #     delta_weddy = weddy - load_phys["weddy"]
    #     delta_wgen = wgen - load_phys["wgen"]
    #     delta_wind = wind - load_phys["wind"]

    #     x = 1

    # if iter == 3719:
    #     conc_name = f"conc_iter{iter:04d}.npy"
    #     rates_name = f"rates_iter{iter:04d}.npy"
    #     phys_name = f"phys_iter{iter:04d}.npz"

    #     load_conc = np.load(os.getcwd() + "/tests/bfm17/check_conc/" + conc_name, allow_pickle=True)
    #     load_rates = np.load(os.getcwd() + "/tests/bfm17/check_rates/" + rates_name, allow_pickle=True)
    #     load_phys = np.load(os.getcwd() + "/tests/bfm17/check_phys/" + phys_name, allow_pickle=True)

    #     delta_conc = conc - load_conc
    #     delta_rates = d_dt - load_rates
    #     delta_temp = temp - load_phys["temp"]
    #     delta_sal = sal - load_phys["sal"]
    #     delta_dens = dens - load_phys["dens"]
    #     delta_ism = ism - load_phys["ism"]
    #     delta_z = z - load_phys["z"]
    #     delta_dz = dz - load_phys["dz"]
    #     delta_surface_PAR = surface_PAR - load_phys["surface_PAR"]
    #     delta_weddy = weddy - load_phys["weddy"]
    #     delta_wgen = wgen - load_phys["wgen"]
    #     delta_wind = wind - load_phys["wind"]

    #     x = 1


    # concentration = vertical_diffusivity(iter, concentration, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
    #                      nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
    #                      vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, umolbgc, nbcbgc, ntp, shortwave_radiation, dif_trac)

    conc_bwd, conc_cur = vertical_diffusivity(iter, conc_bwd, conc_cur, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
                            nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
                            vertical_grid, vertical_spacing, vertical_spacing_staggered, vertical_spacing_reciprocal, umolbgc, nbcbgc, ntp, shortwave_radiation, dif_trac)

    return conc_bwd, conc_cur


@njit
def bgc_physical(temperature, salinity, density, inorganic_suspended_matter, vertical_spacing, column_depth, swrad, water_specific_heat_times_density, w_eddy_velocity, w_gen, wsu, wsv):

    # phys_ids = List.empty_list(unicode_type)
    # phys_params = List.empty_list(float64[:])

    temp = temperature[:-1]
    sal = salinity[:-1]
    dens = (density[:-1] * 1.E+03) + 1.E+03
    ism = inorganic_suspended_matter[:]
    z = vertical_spacing[:-1] * column_depth
    dz = vertical_spacing[:-1]
    surface_PAR = -swrad * water_specific_heat_times_density
    weddy = w_eddy_velocity[:]
    wgen = w_gen[:]

    rms_wind = np.sqrt(wsu**2 + wsv**2) * 1.E+03
    wind = np.sqrt(rms_wind/(1.25 * 0.0014))

    return temp, sal, dens, ism, z, dz, surface_PAR, weddy, wgen, wind


@njit
def vertical_advection(b_cur, b_bwd, b_fwd, sinking_velocity, num_layers, dzr):
    """"
    Description: Handles the sinking of BFM state variablles. Sinking is treated as downward vertical advection
                 computed with upstream finite differences.
    NOTE: Downward velocities are negative
    """
    # sinking velocity input from vdiff_SOS
    b_cur[-1] = b_cur[-2]
    b_bwd[-1] = b_bwd[-2]
    
    b_fwd[0] = dzr[0] * b_cur[0] * sinking_velocity[1]
    for i in range(1,num_layers-1):
        b_fwd[i] = dzr[i] * (b_cur[i] * sinking_velocity[i + 1] - b_cur[i - 1] * sinking_velocity[i])

    return b_fwd


# @njit
# def vertical_diffusivity(iter, concentration, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
#                          nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
#                          z, dz, dzz, dzr, umol, nbc, ntp, swrad, kh):

@njit
def vertical_diffusivity(iter, conc_bwd, conc_cur, d_dt, dt2, num_layers, column_depth, smoth, sinking, weddy, wgen, tracer_map, tracer_type,
                         nrt_o2, nrt_po4, nrt_no3, nrt_nh4, d_o2surf, o2b, no3b, ponb_grad, po4b,
                         z, dz, dzz, dzr, umol, nbc, ntp, swrad, kh):
    """
    Description: Calculates the vertical diffusivity of BFM biochemical components and
                 integrates BFM state variables with Source Splitting (SoS) method
    """
    # Reverse the tracer map
    reverse_map = reverse_tracer_map(tracer_map)

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0

    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/86400.  # to m/s

    # Relaxation velocities
    trelax_o2 = nrt_o2 / 86400.
    trelax_po4 = nrt_po4 / 86400.
    trelax_no3 = nrt_no3 / 86400.
    trelax_nh4 = nrt_nh4

    # if iter > 0:
    #     dif_cur_pre = conc_cur - check_after_vdiff_SOS.conc_cur[iter-1]
    #     dif_bwd_pre = conc_bwd - check_after_vdiff_SOS.conc_bwd[iter-1]


    # if iter >=24 and iter <=29: # matches up to iter==27, wrong starting at iter==28
    #     # concentrations are correct at coming out of iter==27
    #     # rate calculations are incorrect starting at iter==28, cascades through simulation from there
    #     # what happens at iter==28?
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter24to29[iter-24]
    #     x=1
    # x = 1
    # Loop over bgc state variables
    # for i in range(0, len(concentration)):   # i = tracer constituent
    for i in range(0, len(conc_cur)):   # i = tracer constituent
        # Zeroing of previous tracer
        b_cur = np.zeros(num_layers, dtype=np.float64)
        b_bwd = np.zeros(num_layers, dtype=np.float64)
        b_fwd = np.zeros(num_layers, dtype=np.float64)
        b_surf = 0.
        b_sflx = 0.
        b_bflx = 0.

        # Load BFM state variable
        # if iter == 0:   # Initialize backward time level on first iteration
        #     b_cur[:-1] = concentration[i,:,iter]
        #     b_bwd[:-1] = concentration[i,:,iter]
        # else:   # Current and backward time levels previously calculated
        #     b_cur[:-1] = concentration[i,:,iter]
        #     b_bwd[:-1] = concentration[i,:,iter-1]
        b_cur[:-1] = conc_cur[i]
        b_bwd[:-1] = conc_bwd[i]

        b_cur[-1] = b_cur[-2]
        b_bwd[-1] = b_bwd[-2]

        # Calculate tracer sinking velocity
        sinking_velocity = W_ON*wgen + Weddy_ON*weddy

        # if reverse_map[i] == 'o2':
        #     b_sflx = -(d_o2surf / 86400.)
        #     b_bflx = (concentration[i,-1,iter] - o2b) * trelax_o2
        # elif reverse_map[i] == 'no3':
        #     b_sflx = 0.
        #     b_bflx = (concentration[i,-1,iter] - no3b) * trelax_no3
        # elif reverse_map[i] == 'nh4':
        #     b_sflx = 0.
        #     b_bflx = ponb_grad * trelax_nh4
        # elif reverse_map[i] == 'po4':
        #     b_sflx = 0.
        #     b_bflx = (concentration[i,-1,iter] - po4b) * trelax_po4
        # elif reverse_map[i] == 'co2':
        #     b_sflx = 0.
        # elif reverse_map[i] == 'sio4':
        #     b_sflx = 0.
        if reverse_map[i] == 'o2':
            b_sflx = -(d_o2surf / 86400.)
            b_bflx = (conc_cur[i,-1] - o2b) * trelax_o2
        elif reverse_map[i] == 'no3':
            b_sflx = 0.
            b_bflx = (conc_cur[i,-1] - no3b) * trelax_no3
        elif reverse_map[i] == 'nh4':
            b_sflx = 0.
            b_bflx = ponb_grad * trelax_nh4
        elif reverse_map[i] == 'po4':
            b_sflx = 0.
            b_bflx = (conc_cur[i,-1] - po4b) * trelax_po4
        elif reverse_map[i] == 'co2':
            b_sflx = 0.
        elif reverse_map[i] == 'sio4':
            b_sflx = 0.
        
        sinking_velocity[:-1] -= sinking[i] / 86400.
        if tracer_type[i] == "phytoplankton":  
            # Final sink value for phytoplankton
            sinking_velocity[-1] = sinking_velocity[-2]

        if tracer_type[i] == "particulate":
            # Final sink value for particulate detritus
            sinking_velocity[-1] = sinking_velocity[-2]

        # Sinking: upstream vertical advection
        b_fwd = vertical_advection(b_cur, b_bwd, b_fwd, sinking_velocity, num_layers, dzr)
        
        # Source splitting (SoS) leapfrog integration
        for j in range(0,num_layers-1):
            b_fwd[j] = b_bwd[j] + dt2*( (b_fwd[j]/column_depth) + d_dt[i,j] ) #+ tracers[key].d_dt[index,i])
        
        # Compute vertical diffusion and terminate integration (implicit leapfrogging)
        b_fwd, b_surf, b_sflx, b_bflx = temperature_and_salinity_profiles('BGC', dt2, num_layers, column_depth, z, dz, dzz, umol, nbc, ntp, swrad, kh, b_fwd, b_surf, b_sflx, b_bflx)
        
        # Clipping (if needed)
        for j in range(0,num_layers-1):
            b_fwd[j] = max(1.E-20,b_fwd[j])
        
        # Mix the time step and restore time sequence
        # concentration[i,:,iter] = b_cur[:-1] + 0.5 * smoth * (b_fwd[:-1] + b_bwd[:-1] - 2.*b_cur[:-1])
        # concentration[i,:,iter+1] = b_fwd[:-1]

        conc_bwd[i,:] = b_cur[:-1] + 0.5 * smoth * (b_fwd[:-1] + b_bwd[:-1] - 2.*b_cur[:-1])
        conc_cur[i,:] = b_fwd[:-1]

    # if iter < 5:
    #     # incorrect rates at iter==0 for
    #     # bac1 [c,n,p] [7,8,9]
    #     # dom1 [c,n,p] [39,40,41]
    #     # pom1 [c,n,p] [44,45,46] (pom1 [s] is fine)
    #     # rates are lower than expected (dif_ddt = +) for bac1, higher than expected  (dif_ddt = -) for dom1 and pom1
    #     # all other rates being correct likely points to bacteria uptake or mortality
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur[iter]
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd[iter]
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates[iter]
    #     x=1

    # if iter == 24:  # matches here
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter24
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter24
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter24
    #     x=1

    # if iter >=24 and iter <=29: # matches up to iter==27, wrong starting at iter==28
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter24to29[iter-24]
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter24to29[iter-24]
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter24to29[iter-24]
    #     x=1

    # if iter == 28:
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter28
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter28
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter28
    #     x=1

    # if iter == 29:
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter29
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter29
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter29
    #     x=1

    # if iter == 49:  # wrong by here, fixed (i think) by turning on sedimentation for phytoplankton
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter49
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter49
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter49
    #     x=1

    # if iter == 119:
    #     dif_cur = conc_cur - check_after_vdiff_SOS.conc_cur_iter119
    #     dif_bwd = conc_bwd - check_after_vdiff_SOS.conc_bwd_iter119
    #     dif_ddt = d_dt - check_after_vdiff_SOS.rates_iter119
    #     x=1
    
    # return concentration
    return conc_bwd, conc_cur


@njit
def get_tracer_from_index(index, tracer_map):
    """
    Definition: Identifies the tracer key associated with index
    :return: tracer key
    """
    for key,value in tracer_map.items():
        if index in value:
            tracer = key

    return tracer


@njit
def reverse_tracer_map(tracer_map):
    """
    Definition: Reverses tracer map for quicker index lookup
    :return: reversed tracer map
    """
    reverse_map = Dict.empty(key_type=types.int64, value_type=types.unicode_type)

    for tracer,value in tracer_map.items():
        for index in value:
            reverse_map[index] = tracer

    return reverse_map
