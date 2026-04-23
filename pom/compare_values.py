from pom.check_rates import o2rates, no3rates, nh4rates, po4rates, phytocrates, phytonrates, phytoprates, phytolrates, zoocrates, zoonrates, zooprates, domcrates, domnrates, domprates, pomcrates, pomnrates, pomprates
from pom.check_conc import o2, no3, nh4, po4, phyto1, zoo1, dom1, pom1

def compare_values(iter, tracers):

    if iter < 3:
        d_dt_diff = {
            "o2": tracers["o2"].d_dt[0,:] - o2rates[iter],
            "no3": tracers["no3"].d_dt[0,:] - no3rates[iter],
            "nh4": tracers["nh4"].d_dt[0,:] - nh4rates[iter],
            "po4": tracers["po4"].d_dt[0,:] - po4rates[iter],
            "phyc": tracers["phyto1"].d_dt[0,:] - phytocrates[iter],
            "phyn": tracers["phyto1"].d_dt[1,:] - phytonrates[iter],
            "phyp": tracers["phyto1"].d_dt[2,:] - phytoprates[iter],
            "phyl": tracers["phyto1"].d_dt[3,:] - phytolrates[iter],
            "zooc": tracers["zoo1"].d_dt[0,:] - zoocrates[iter],
            "zoon": tracers["zoo1"].d_dt[1,:] - zoonrates[iter],
            "zoop": tracers["zoo1"].d_dt[2,:] - zooprates[iter],
            "domc": tracers["dom1"].d_dt[0,:] - domcrates[iter],
            "domn": tracers["dom1"].d_dt[1,:] - domnrates[iter],
            "domp": tracers["dom1"].d_dt[2,:] - domprates[iter],
            "pomc": tracers["pom1"].d_dt[0,:] - pomcrates[iter],
            "pomn": tracers["pom1"].d_dt[1,:] - pomnrates[iter],
            "pomp": tracers["pom1"].d_dt[2,:] - pomprates[iter]
        }
    else:
        d_dt_diff = {
            "o2": None,
            "no3": None,
            "nh4": None,
            "po4": None,
            "phyc": None,
            "phyn": None,
            "phyp": None,
            "phyl": None,
            "zooc": None,
            "zoon": None,
            "zoop": None,
            "domc": None,
            "domn": None,
            "domp": None,
            "pomc": None,
            "pomn": None,
            "pomp": None
        }

    conc_diff = {
        "o2": tracers["o2"].conc[0,:,iter] - o2[iter],
        "no3": tracers["no3"].conc[0,:,iter] - no3[iter],
        "nh4": tracers["nh4"].conc[0,:,iter] - nh4[iter],
        "po4": tracers["po4"].conc[0,:,iter] - po4[iter],
        "phyc": tracers["phyto1"].conc[0,:,iter] - phyto1[0][iter],
        "phyn": tracers["phyto1"].conc[1,:,iter] - phyto1[1][iter],
        "phyp": tracers["phyto1"].conc[2,:,iter] - phyto1[2][iter],
        "phyl": tracers["phyto1"].conc[3,:,iter] - phyto1[3][iter],
        "zooc": tracers["zoo1"].conc[0,:,iter] - zoo1[0][iter],
        "zoon": tracers["zoo1"].conc[1,:,iter] - zoo1[1][iter],
        "zoop": tracers["zoo1"].conc[2,:,iter] - zoo1[2][iter],
        "domc": tracers["dom1"].conc[0,:,iter] - dom1[0][iter],
        "domn": tracers["dom1"].conc[1,:,iter] - dom1[1][iter],
        "domp": tracers["dom1"].conc[2,:,iter] - dom1[2][iter],
        "pomc": tracers["pom1"].conc[0,:,iter] - pom1[0][iter],
        "pomn": tracers["pom1"].conc[1,:,iter] - pom1[1][iter],
        "pomp": tracers["pom1"].conc[2,:,iter] - pom1[2][iter],
    }

    return d_dt_diff, conc_diff