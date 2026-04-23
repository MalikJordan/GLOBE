import numpy as np
from pom.calculations import temperature_and_salinity_profiles
from functions.bgc_rate_eqns import bgc_rate_eqns
import os
# from pom.check_rates import o2rates, no3rates, nh4rates, po4rates, phytocrates, phytonrates, phytoprates, phytolrates, zoocrates, zoonrates, zooprates, domcrates, domnrates, domprates, pomcrates, pomnrates, pomprates
# from pom.check_conc import o2, no3, nh4, po4, phyto1, zoo1, dom1, pom1
from pom.compare_values import compare_values

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

def pom_bgc_1d(iter, base_element, physical, pom1d, tracers):

    physical = bgc_physical(physical, pom1d)
    bgc_rate_eqns(iter, base_element, physical, pom1d, tracers)
    # print('i= ',iter)
    # # print(tracers["phyto1"].conc[0,0,iter] - phytoc0[iter])
    # # print(tracers["phyto1"].d_dt[0,0] - dphytoc0_dt[iter])

    # globe_rates[iter,0] = tracers["o2"].d_dt[0,0]
    # globe_rates[iter,1] = tracers["po4"].d_dt[0,0]
    # globe_rates[iter,2] = tracers["no3"].d_dt[0,0]
    # globe_rates[iter,3] = tracers["nh4"].d_dt[0,0]
    # globe_rates[iter,4] = tracers["hs"].d_dt[0,0]
    # globe_rates[iter,5] = tracers["phyto1"].d_dt[0,0]
    # globe_rates[iter,6] = tracers["phyto1"].d_dt[1,0]
    # globe_rates[iter,7] = tracers["phyto1"].d_dt[2,0]
    # globe_rates[iter,8] = tracers["phyto1"].d_dt[3,0]
    # globe_rates[iter,9] = tracers["zoo1"].d_dt[0,0]
    # globe_rates[iter,10] = tracers["zoo1"].d_dt[1,0]
    # globe_rates[iter,11] = tracers["zoo1"].d_dt[2,0]
    # globe_rates[iter,12] = tracers["dom1"].d_dt[0,0]
    # globe_rates[iter,13] = tracers["dom1"].d_dt[1,0]
    # globe_rates[iter,14] = tracers["dom1"].d_dt[2,0]
    # globe_rates[iter,15] = tracers["pom1"].d_dt[0,0]
    # globe_rates[iter,16] = tracers["pom1"].d_dt[1,0]
    # globe_rates[iter,17] = tracers["pom1"].d_dt[2,0]

    # print(globe_rates[iter,:] - bfm_rates[iter])
    o2_surf = tracers["o2"].conc[0,0,iter]
    no3_surf = tracers["no3"].conc[0,0,iter]
    nh4_surf= tracers["nh4"].conc[0,0,iter]
    po4_surf = tracers["po4"].conc[0,0,iter]
    phyc_surf = tracers["phyto1"].conc[0,0,iter]
    phyn_surf = tracers["phyto1"].conc[1,0,iter]
    phyp_surf = tracers["phyto1"].conc[2,0,iter]
    phyl_surf = tracers["phyto1"].conc[3,0,iter]
    zc_surf = tracers["zoo1"].conc[0,0,iter]
    zn_surf = tracers["zoo1"].conc[1,0,iter]
    zp_surf = tracers["zoo1"].conc[2,0,iter]
    dc_surf = tracers["dom1"].conc[0,0,iter]
    dn_surf = tracers["dom1"].conc[1,0,iter]
    dp_surf = tracers["dom1"].conc[2,0,iter]
    pc_surf = tracers["pom1"].conc[0,0,iter]
    pn_surf = tracers["pom1"].conc[1,0,iter]
    pp_surf = tracers["pom1"].conc[2,0,iter]
    o2_bot = tracers["o2"].conc[0,-1,iter]
    no3_bot = tracers["no3"].conc[0,-1,iter]
    nh4_bot = tracers["nh4"].conc[0,-1,iter]
    po4_bot = tracers["po4"].conc[0,-1,iter]
    phyc_bot = tracers["phyto1"].conc[0,-1,iter]
    phyn_bot = tracers["phyto1"].conc[1,-1,iter]
    phyp_bot = tracers["phyto1"].conc[2,-1,iter]
    phyl_bot = tracers["phyto1"].conc[3,-1,iter]
    zc_bot = tracers["zoo1"].conc[0,-1,iter]
    zn_bot = tracers["zoo1"].conc[1,-1,iter]
    zp_bot = tracers["zoo1"].conc[2,-1,iter]
    dc_bot = tracers["dom1"].conc[0,-1,iter]
    dn_bot = tracers["dom1"].conc[1,-1,iter]
    dp_bot = tracers["dom1"].conc[2,-1,iter]
    pc_bot = tracers["pom1"].conc[0,-1,iter]
    pn_bot = tracers["pom1"].conc[1,-1,iter]
    pp_bot = tracers["pom1"].conc[2,-1,iter]

    if iter < 5:
        d_dt_diff, conc_diff = compare_values(iter, tracers)
    # if iter == 9:
    #     x=1
    vertical_diffusivity(iter, physical, pom1d, tracers)


def bgc_physical(physical, pom1d):

    # Initialize dictionary
    physical["bgc_phys_vars"] = {}

    # 1D arrays for bgc calculations
    physical["bgc_phys_vars"]["temperature"] = physical["temperature"]["tb"][:-1]
    physical["bgc_phys_vars"]["salinity"] = physical["salinity"]["sb"][:-1]
    physical["bgc_phys_vars"]["density"] = (physical["density"][:-1] * 1.E+03) + 1.E+03
    physical["bgc_phys_vars"]["ism"] = physical["ism"]
    physical["bgc_phys_vars"]["z"] = physical["vertical_grid"]["dz"][:-1] * physical["water_column"]["column_depth"]
    physical["bgc_phys_vars"]["dz"] = physical["vertical_grid"]["dz"][:-1]
    physical["bgc_phys_vars"]["surface_PAR"] = -physical["swrad"] * pom1d["general"]["water_specific_heat_times_density"]
    physical["bgc_phys_vars"]["weddy"] = physical["weddy"]
    physical["bgc_phys_vars"]["wgen"] = physical["wgen"]

    wind = np.sqrt(physical["stresses"]["wsu"]**2 + physical["stresses"]["wsv"]**2) * 1.E+03
    physical["bgc_phys_vars"]["wind"] = np.sqrt(wind/(1.25 * 0.0014))

    return physical


def vertical_advection(physical, bgc_state_var, sinking_velocity):
    """"
    Description: Handles the sinking of BFM state variablles. Sinking is treated as downward vertical advection
                 computed with upstream finite differences.
    NOTE: Downward velocities are negative
    """
    # sinking velocity input from vdiff_SOS
    bgc_state_var["b"][-1] = bgc_state_var["b"][-2]
    bgc_state_var["bb"][-1] = bgc_state_var["bb"][-2]
    
    bgc_state_var["bf"][0] = physical["vertical_grid"]["dzr"][0] * bgc_state_var["b"][0] * sinking_velocity[1]
    for i in range(1,physical["water_column"]["num_layers"]-1):
        bgc_state_var["bf"][i] = physical["vertical_grid"]["dzr"][i] * (bgc_state_var["b"][i] * sinking_velocity[i + 1] - bgc_state_var["b"][i - 1] * sinking_velocity[i])

    return bgc_state_var


def vertical_diffusivity(iter, physical, pom1d, tracers):
    """
    Description: Calculates the vertical diffusivity of BFM biochemical components and
                 integrats BFM state variables with Source Splitting (SoS) method
    """

    # The input general cir. vertical vel. is suppose to be in m/s
    W_ON = 1.0

    # The input eddy vertical vel. is provided in m/d
    Weddy_ON = 0.1/physical["simulation"]["sec_per_day"]  # to m/s

    # Relaxation velocities
    trelax_o2 = pom1d["relaxation_velocities"]["nrt_o2"] / physical["simulation"]["sec_per_day"]
    trelax_po4 = pom1d["relaxation_velocities"]["nrt_po4"] / physical["simulation"]["sec_per_day"]
    trelax_no3 = pom1d["relaxation_velocities"]["nrt_no3"] / physical["simulation"]["sec_per_day"]
    trelax_nh4 = pom1d["relaxation_velocities"]["nrt_nh4"]

    # Loop over bgc state variables
    for key in tracers:
        if key == "phyto1":
            x=1
        for const in tracers[key].composition:
            # Get tracer for time stepping
            index = tracers[key].composition.index(const)
            
            # Zeroing of previous tracer
            bgc_state_var = {
                "b": np.zeros(physical["water_column"]["num_layers"]),  # Current
                "bf": np.zeros(physical["water_column"]["num_layers"]), # Forward
                "bb": np.zeros(physical["water_column"]["num_layers"]), # Backward
                "surf": 0.,         # Surface value
                "surf_flux": 0.,    # Surface flux
                "bot_flux": 0.      # Bottom flux
            }

            # Load BFM state variable
            if iter == 0:   # Initialize backward time level on first iteration
                bgc_state_var["b"][:-1] = tracers[key].conc[index,:,iter]
                bgc_state_var["bb"][:-1] = tracers[key].conc[index,:,iter]
            else:   # Current and backward time levels previously calculated
                bgc_state_var["b"][:-1] = tracers[key].conc[index,:,iter]
                bgc_state_var["bb"][:-1] = tracers[key].conc[index,:,iter-1]
            
            bgc_state_var["b"][-1] = bgc_state_var["b"][-2]
            bgc_state_var["bb"][-1] = bgc_state_var["bb"][-2]

            # Calculate tracer sinking velocity
            sinking_velocity = W_ON*physical["bgc_phys_vars"]["wgen"] + Weddy_ON*physical["bgc_phys_vars"]["weddy"]
            
            if key == 'o2':
                bgc_state_var["surf_flux"] = -(tracers[key].surf_flux[0,0] / physical["simulation"]["sec_per_day"])
                bgc_state_var["bot_flux"] = (tracers[key].conc[index,-1,iter] - physical["nutrients"]["o2b"]) * trelax_o2
            elif key == 'no3':
                bgc_state_var["surf_flux"] = 0.
                bgc_state_var["bot_flux"] = (tracers[key].conc[index,-1,iter] - physical["nutrients"]["no3b"]) * trelax_no3
            elif key == 'nh4':
                bgc_state_var["surf_flux"] = 0.
                bgc_state_var["bot_flux"] = physical["nutrients"]["ponb_grad"] * trelax_nh4
            elif key == 'po4':
                bgc_state_var["surf_flux"] = 0.
                bgc_state_var["bot_flux"] = (tracers[key].conc[index,-1,iter] - physical["nutrients"]["po4b"]) * trelax_po4
            # elif key == 'co2':
            #     bgc_state_var["surf_flux"] = 0.
            #     bgc_state_var["bot_flux"] = (tracers[key].conc[index,-1,iter] - physical["nutrients"]["no3b"]) * trelax_no3
            # elif key == 'sio4':
            #     bgc_state_var["surf_flux"] = 0.
            #     bgc_state_var["bot_flux"] = (tracers[key].conc[index,-1,iter] - physical["nutrients"]["no3b"]) * trelax_no3

            if hasattr(tracers[key],"sinking_velocity"):
                # Include additional sinking velocity for sinking tracers
                sinking_velocity[:-1] -= tracers[key].sinking_velocity / physical["simulation"]["sec_per_day"]
                
            if tracers[key].type == "phytoplankton":  
                # Final sink value for phytoplankton
                sinking_velocity[-1] = sinking_velocity[-2]

            if tracers[key].type == "detritus":
                # Final sink value for particulate detritus
                if tracers[key].form == "particulate":  sinking_velocity[-1] = sinking_velocity[-2]

            # Sinking: upstream vertical advection
            bgc_state_var = vertical_advection(physical, bgc_state_var, sinking_velocity)
            
            # Source splitting (SoS) leapfrog integration
            for i in range(0,physical["water_column"]["num_layers"]-1):
                bgc_state_var["bf"][i] = bgc_state_var["bb"][i] + physical["simulation"]["dt2"]*((bgc_state_var["bf"][i]/physical["water_column"]["column_depth"]) + tracers[key].d_dt[index,i])
            
            # Compute vertical diffusion and terminate integration (implicit leapfrogging)
            bgc_state_var = temperature_and_salinity_profiles(physical, pom1d, bgc_state_var, 'BGC')
            
            if key == 'o2':
                x=1
            
            # Clipping (if needed)
            for i in range(0,physical["water_column"]["num_layers"]-1):
                bgc_state_var["bf"][i] = max(1.E-20,bgc_state_var["bf"][i])
            
            # Mix the time step and restore time sequence
            tracers[key].conc[index,:,iter] = bgc_state_var["b"][:-1] + 0.5 * pom1d["general"]["smoth"] * (bgc_state_var["bf"][:-1] + bgc_state_var["bb"][:-1] - 2.*bgc_state_var["b"][:-1])
            tracers[key].conc[index,:,iter+1] = bgc_state_var["bf"][:-1]

    return 

