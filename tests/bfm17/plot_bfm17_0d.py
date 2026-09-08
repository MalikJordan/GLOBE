import netCDF4 as nc
import os
import numpy as np
from matplotlib import pyplot as plt

folder = os.getcwd() + '/tests/bfm17'

# Get Fortran Data --------------------------------------------------------------------
path = os.getcwd() + '/tests/bfm17/data/BFM_standalone_pelagic.nc'
variables = nc.Dataset(path).variables

o2_bfm17 = np.asarray(variables['O2o'][:360])
po4_bfm17 = np.asarray(variables['N1p'][:360])
no3_bfm17 = np.asarray(variables['N3n'][:360])
nh4_bfm17 = np.asarray(variables['N4n'][:360])

pc_bfm17 = np.asarray(variables['P2c'][:360])
domc_bfm17 = np.asarray(variables['R1c'][:360])
pomc_bfm17 = np.asarray(variables['R6c'][:360])
zc_bfm17 = np.asarray(variables['Z5c'][:360])

# Get GLOBE Data --------------------------------------------------------------------
# Load solution
# path = os.getcwd() + "/tests/bfm17/data/bfm17.npz"
# path = os.getcwd() + "/concentration_bfm17_0d.npz"
path = os.getcwd() + "/concentration_bfm17_0d_0907.npz"
model = np.load(path, allow_pickle=True)
globe_daily = model["daily"]
globe_monthly = model["monthly"]

globe_daily = globe_daily[...,:360]
globe_monthly = globe_monthly[...,:12]

# Load tracer indices
# path = os.getcwd() + "/tests/bfm17/data/tracer_indices_bfm17.npz"
path = os.getcwd() + "/tracer_indices_bfm17_0d.npz"
path = os.getcwd() + "/tracer_indices_bfm17_0d_0907.npz"
indices = np.load(path)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

o2 = tracer_indices["o2"][0]
no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
po4 = tracer_indices["po4"][0]
pc = tracer_indices["phyto1"][0]
pn = tracer_indices["phyto1"][1]
pp = tracer_indices["phyto1"][2]
pl = tracer_indices["phyto1"][3]
zc = tracer_indices["zoo1"][0]
zn = tracer_indices["zoo1"][1]
zp = tracer_indices["zoo1"][2]
domc = tracer_indices["dom1"][0]
domn = tracer_indices["dom1"][1]
domp = tracer_indices["dom1"][2]
pomc = tracer_indices["pom1"][0]
pomn = tracer_indices["pom1"][1]
pomp = tracer_indices["pom1"][2]

# xlabel = ['J','A','J','O','J','A','J','O','J']
# xticks = [0,90,180,270,360,450,540,630,720]
xlabel = ['J','A','J','O','J']
xticks = [0,90,180,270,360]
days = np.linspace(0,359,360)

# Create Plots --------------------------------------------------------------------
# Nutrients
fig, axs = plt.subplots(2,2,figsize=(10,10),sharex=True)

axs[0,0].plot(days,globe_daily[o2,0],label='GLOBE')
axs[0,0].plot(days,o2_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[0,0].set_title("(a) Oxygen")
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mmol O ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[no3,0],label='GLOBE')
axs[0,1].plot(days,no3_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[0,1].set_title("(b) Nitrate")
axs[0,1].set_xlim([0,360])
axs[0,1].set_ylabel("mmol N ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[nh4,0],label='GLOBE')
axs[1,0].plot(days,nh4_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[1,0].set_title("(c) Ammonium")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mmol N ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[po4,0],label='GLOBE')
axs[1,1].plot(days,po4_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[1,1].set_title("(d) Phosphate")
axs[1,1].set_xlabel("Time [months]")
axs[1,1].set_xticks(xticks,xlabel)
axs[1,1].set_xlim([0,360])
axs[1,1].set_ylabel("mmol P ${m^{-3}}$")

handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles,labels, loc='lower center', ncol=2)

# fig.suptitle("Nutrients")
fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
nut = os.path.join(folder + "/figures","nutrients.jpg")
plt.savefig(nut)


# Organic
fig, axs = plt.subplots(2,2,figsize=(10,10), sharex=True)

axs[0,0].plot(days,globe_daily[pc,0],label='GLOBE')
axs[0,0].plot(days,pc_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[0,0].set_title("(a) Phytoplankton -- Carbon")
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mg C ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[zc,0],label='GLOBE')
axs[0,1].plot(days,zc_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[0,1].set_title("(b) Zooplankton -- Carbon")
axs[0,1].set_xlim([0,360])
axs[0,1].set_ylabel("mg C ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[domc,0],label='GLOBE')
axs[1,0].plot(days,domc_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[1,0].set_title("(c) Dissolved Organic Carbon")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mg C ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[pomc,0],label='GLOBE')
axs[1,1].plot(days,pomc_bfm17,linestyle=(0, (5, 10)),color='black',label='BFM17')
axs[1,1].set_title("(d) Particulate Organic Carbon")
axs[1,1].set_xlabel("Time [months]")
axs[1,1].set_xticks(xticks,xlabel)
axs[1,1].set_xlim([0,360])
axs[1,1].set_ylabel("mg C ${m^{-3}}$")

handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles,labels, loc='lower center', ncol=2)

# fig.suptitle("Organic Carbon")
fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
org = os.path.join(folder + "/figures","organic.jpg")
plt.savefig(org)