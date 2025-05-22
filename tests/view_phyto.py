import numpy as np
import netCDF4 as nc
from matplotlib import pyplot as plt
import os


def check_gpp():

    path = os.getcwd() + '/bfm17/BFM_standalone_pelagic.nc'
    variables = nc.Dataset(path).variables

    exu_bfm17 = np.asarray(variables['exPPYc'][:])
    gpp_bfm17 = np.asarray(variables['ruPPYc'][:])

    phyto = np.load("phyto-0424.npz", allow_pickle=True)
    gpp = phyto['gpp']
    npp = phyto['npp']
    psn = phyto['psn']
    rsp = phyto['rsp']
    exu = phyto['exu']

    iters = np.linspace(0,len(gpp)-1,len(gpp))

    iters_per_day = 86400/360
    gpp_avg = np.zeros(int(np.floor((len(iters))/iters_per_day)))
    npp_avg = np.zeros_like(gpp_avg)
    psn_avg = np.zeros_like(gpp_avg)
    rsp_avg = np.zeros_like(gpp_avg)
    exu_avg = np.zeros_like(gpp_avg)

    day=0
    for i in range(0,len(iters)):
        gpp_avg[day] += gpp[i]
        npp_avg[day] += npp[i]
        psn_avg[day] += psn[i]
        rsp_avg[day] += rsp[i]
        exu_avg[day] += exu[i]

        if i%int(iters_per_day) == 0 and i != 0:
            gpp_avg[day] = gpp_avg[day]/iters_per_day
            npp_avg[day] = npp_avg[day]/iters_per_day
            psn_avg[day] = psn_avg[day]/iters_per_day
            rsp_avg[day] = rsp_avg[day]/iters_per_day
            exu_avg[day] = exu_avg[day]/iters_per_day

            day += 1
            if day == 730:
                break

    fig,ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),gpp_bfm17[:730],'k')
    ax.plot(np.linspace(0,729,730),gpp_avg,'--r')
    ax.set_title('gross primary production')
    ax.legend(['bfm17','globe'])
    plt.savefig('gpp-2.jpg')

    fig,ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),exu_bfm17[:730],'k')
    ax.plot(np.linspace(0,729,730),exu_avg,'--r')
    ax.set_title('carbon excretion')
    ax.legend(['bfm17','globe'])
    plt.savefig('exu-2.jpg')


def check_nut_lim():

    nut = np.load("nutrient_limitation.npz", allow_pickle=True)
    nut_lim = nut["nut_lim"]
    nut_globe = np.load("nut_lim_globe.npz", allow_pickle=True)
    nut_lim_globe = nut_globe["nut_lim"]

    iters = np.linspace(0,len(nut_lim_globe)-1,len(nut_lim_globe))
    iters_per_day = 86400/360
    nut_avg = np.zeros(int(np.floor((len(iters))/iters_per_day)))
    nut_avg_globe = np.zeros_like(nut_avg)

    day=0
    for i in range(0,len(iters)):
        nut_avg[day] += nut_lim[i]
        nut_avg_globe[day] += nut_lim_globe[i]

        if i%int(iters_per_day) == 0 and i != 0:
            nut_avg[day] = nut_avg[day]/iters_per_day
            nut_avg_globe[day] = nut_avg_globe[day]/iters_per_day

            day += 1
            if day == 730:
                break

    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),nut_avg,'k')
    ax.plot(np.linspace(0,729,730),nut_avg_globe,'--r')
    # ax.legend()
    ax.set_title('nutrient limitation')
    plt.savefig('nut_lim-2.jpg')



def check_light_lim():

    path = os.getcwd() + '/bfm17/BFM_standalone_pelagic.nc'
    variables = nc.Dataset(path).variables

    fI = np.asarray(variables['eiPPY_iiP2_'][:])
    irrad = np.asarray(variables['EIR'][:])

    light_globe = np.load("light_globe.npz", allow_pickle=True)
    irrad_globe = light_globe["irrad"]
    light_lim_globe = light_globe["light_lim"]

    iters = np.linspace(0,len(light_lim_globe)-1,len(light_lim_globe))
    iters_per_day = 86400/360
    irrad_avg = np.zeros(int(np.floor((len(iters))/iters_per_day)))
    irrad_avg_globe = np.zeros_like(irrad_avg)
    light_lim_avg = np.zeros_like(irrad_avg)
    light_lim_avg_globe = np.zeros_like(irrad_avg)

    day=0
    for i in range(0,len(iters)):
        irrad_avg_globe[day] += irrad_globe[i]
        light_lim_avg_globe[day] += light_lim_globe[i]

        if i%int(iters_per_day) == 0 and i != 0:
            irrad_avg_globe[day] = irrad_avg_globe[day]/iters_per_day
            light_lim_avg_globe[day] = light_lim_avg_globe[day]/iters_per_day

            day += 1
            if day == 730:
                break

    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),irrad[:730],'k')
    ax.plot(np.linspace(0,729,730),irrad_avg_globe,'--r')
    # ax.legend()
    ax.set_title('irradiance')
    plt.savefig('irradiance-2.jpg')

    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),fI[:730],'k')
    ax.plot(np.linspace(0,729,730),light_lim_avg_globe,'--r')
    # ax.legend()
    ax.set_title("light limitation factor")
    plt.savefig('light_lim-2.jpg')


def check_uptake():

    path = os.getcwd() + '/bfm17/BFM_standalone_pelagic.nc'
    variables = nc.Dataset(path).variables

    uptn = np.asarray(variables['ruPPYn'][:])
    uptp = np.asarray(variables['ruPPYp'][:])

    phyto = np.load("phyto-0424.npz", allow_pickle=True)
    uptn_globe = phyto['uptn']
    uptp_globe = phyto['uptp']

    iters = np.linspace(0,len(uptn_globe)-1,len(uptn_globe))
    iters_per_day = 86400/360
    uptn_avg_globe = np.zeros(int(np.floor((len(iters))/iters_per_day)))
    uptp_avg_globe = np.zeros_like(uptn_avg_globe)

    day=0
    for i in range(0,len(iters)):
        uptn_avg_globe[day] += uptn_globe[i]
        uptp_avg_globe[day] += uptp_globe[i]

        if i%int(iters_per_day) == 0 and i != 0:
            uptn_avg_globe[day] = uptn_avg_globe[day]/iters_per_day
            uptp_avg_globe[day] = uptp_avg_globe[day]/iters_per_day

            day += 1
            if day == 730:
                break

    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),uptn[:730],'k')
    ax.plot(np.linspace(0,729,730),uptn_avg_globe,'--r')
    # ax.legend()
    ax.set_title('nitrogen uptake')
    plt.savefig('uptn-2.jpg')

    fig, ax = plt.subplots()
    ax.plot(np.linspace(0,729,730),uptp[:730],'k')
    ax.plot(np.linspace(0,729,730),uptp_avg_globe,'--r')
    # ax.legend()
    ax.set_title('phosphorus uptake')
    plt.savefig('uptp-2.jpg')

    # fig, ax = plt.subplots()
    # ax.plot(np.linspace(0,728,729),uptn[:729],'k')
    # ax.plot(np.linspace(0,728,729),uptn_avg_globe[1:],'--r')
    # # ax.legend()
    # ax.set_title('nitrogen uptake')
    # plt.savefig('uptn.jpg')

    # fig, ax = plt.subplots()
    # ax.plot(np.linspace(0,728,729),uptp[:729],'k')
    # ax.plot(np.linspace(0,728,729),uptp_avg_globe[1:],'--r')
    # # ax.legend()
    # ax.set_title('phosphorus uptake')
    # plt.savefig('uptp.jpg')



check_light_lim()
check_nut_lim()
check_gpp()
check_uptake()