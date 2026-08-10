from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os


def nrmse(check,comp):

    avg = np.zeros(len(check))
    dif = np.zeros(len(check))
    rms = np.zeros(len(check))
    max = np.zeros(len(check))
    for i in range(0,len(check)):
        avg[i] = np.abs(np.mean(check[i,:,:]))
        dif[i] = np.max(check[i,:,:]) - np.min(check[i,:,:])
        rms[i] = np.power( np.mean( np.power( check[i,:,:]-comp[i,:,:], 2 ) )   ,0.5)
        max[i] = np.max(check[i,:,:])
    # nrmse = 100*rms/avg    
    nrmse = 100*rms/max

    return rms, nrmse


# ----------------------------------------------------------------------------------------------------
# Extract model results
# ----------------------------------------------------------------------------------------------------
# GLOBE
# path = os.getcwd() + '/concentration_bfm56-180d.npz'
# path = os.getcwd() + '/concentration_bfm56-2yr.npz'
path = os.getcwd() + '/concentration_bfm56-5yr.npz'
globe_bfm = np.load(path, allow_pickle=True)

globe_daily = globe_bfm["daily"]
globe_monthly = globe_bfm["monthly"]

# globe_daily = globe_daily[:,:,:180]
# globe_monthly = globe_monthly[:,:,:6]
# globe_daily = globe_daily[:,:,:720]
# globe_monthly = globe_monthly[:,:,:24]
globe_daily = globe_daily[:,:,:1800]
globe_monthly = globe_monthly[:,:,:60]

# Tracer Indices
# path = os.getcwd() + '/tracer_indices_bfm56-180d.npz'
# path = os.getcwd() + '/tracer_indices_bfm56-2yr.npz'
path = os.getcwd() + '/tracer_indices_bfm56-5yr.npz'
indices = np.load(path)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

o2 = tracer_indices["o2"][0]
co2 = tracer_indices["co2"][0]
no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
n2 = tracer_indices["n2"][0]
po4 = tracer_indices["po4"][0]
sio4 = tracer_indices["sio4"][0]
hs = tracer_indices["hs"][0]
ta = tracer_indices["ta"][0]

bac1c = tracer_indices["bac1"][0]
bac1n = tracer_indices["bac1"][1]
bac1p = tracer_indices["bac1"][2]

phyto1c = tracer_indices["phyto1"][0]
phyto1n = tracer_indices["phyto1"][1]
phyto1p = tracer_indices["phyto1"][2]
phyto1l = tracer_indices["phyto1"][3]
phyto1s = tracer_indices["phyto1"][4]

phyto2c = tracer_indices["phyto2"][0]
phyto2n = tracer_indices["phyto2"][1]
phyto2p = tracer_indices["phyto2"][2]
phyto2l = tracer_indices["phyto2"][3]

phyto3c = tracer_indices["phyto3"][0]
phyto3n = tracer_indices["phyto3"][1]
phyto3p = tracer_indices["phyto3"][2]
phyto3l = tracer_indices["phyto3"][3]

phyto4c = tracer_indices["phyto4"][0]
phyto4n = tracer_indices["phyto4"][1]
phyto4p = tracer_indices["phyto4"][2]
phyto4l = tracer_indices["phyto4"][3]

mesoz1c = tracer_indices["mesozoo1"][0]
mesoz1n = tracer_indices["mesozoo1"][1]
mesoz1p = tracer_indices["mesozoo1"][2]

mesoz2c = tracer_indices["mesozoo2"][0]
mesoz2n = tracer_indices["mesozoo2"][1]
mesoz2p = tracer_indices["mesozoo2"][2]

microz1c = tracer_indices["microzoo1"][0]
microz1n = tracer_indices["microzoo1"][1]
microz1p = tracer_indices["microzoo1"][2]

microz2c = tracer_indices["microzoo2"][0]
microz2n = tracer_indices["microzoo2"][1]
microz2p = tracer_indices["microzoo2"][2]

dom1c = tracer_indices["dom1"][0]
dom1n = tracer_indices["dom1"][1]
dom1p = tracer_indices["dom1"][2]
dom2c = tracer_indices["dom2"][0]
dom3c = tracer_indices["dom3"][0]

pom1c = tracer_indices["pom1"][0]
pom1n = tracer_indices["pom1"][1]
pom1p = tracer_indices["pom1"][2]
pom1s = tracer_indices["pom1"][3]

# pyPOM1D-BFM
# path = os.getcwd() + '/tests/bfm56/data/pyPOM1D-BFM50-180d.npz'
# path = os.getcwd() + '/tests/bfm56/data/pyPOM1D-BFM50-2yr.npz'
path = os.getcwd() + '/tests/bfm56/data/pyPOM1D-BFM50-5yr.npz'
pypom_bfm = np.load(path, allow_pickle=True)

pypom_daily = pypom_bfm["daily"]
pypom_monthly = pypom_bfm["monthly"]

# pypom_daily = pypom_daily[:,:,:180]
# pypom_monthly = pypom_monthly[:,:,:6]
# pypom_daily = pypom_daily[:,:,:720]
# pypom_monthly = pypom_monthly[:,:,:24]
pypom_daily = pypom_daily[:,:,:1800]
pypom_monthly = pypom_monthly[:,:,:60]


# ----------------------------------------------------------------------------------------------------
# Line Plots
# ----------------------------------------------------------------------------------------------------
# days = np.linspace(0,179,180)
# xticks = [15,45,75,105,135,165]
# xlabel = ['J','','M','','M','']
# days = np.linspace(0,719,720)
# xticks = [90,180,270,360,450,540,630,720]
days = np.linspace(0,1799,1800)
xticks = [180,360,540,720,900,1080,1260,1440,1620,1800]

# Oxygen -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Oxygen [mmol $O_{2}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[o2,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[o2,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[o2,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[o2,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[o2,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[o2,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[o2,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[o2,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[o2,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[o2,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[o2,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[o2,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[o2,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[o2,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/oxygen.jpg')
plt.close()

# Dissolved Inorganic Carbon -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dissolved Inorganic Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[co2,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[co2,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[co2,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[co2,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[co2,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[co2,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[co2,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[co2,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[co2,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[co2,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[co2,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[co2,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[co2,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[co2,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/co2.jpg')
plt.close()

# Nitrate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Nitrate [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[no3,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[no3,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[no3,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[no3,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[no3,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[no3,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[no3,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[no3,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[no3,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[no3,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[no3,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[no3,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[no3,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[no3,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/nitrate.jpg')
plt.close()

# Ammonium -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Ammonium [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[nh4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[nh4,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[nh4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[nh4,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[nh4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[nh4,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[nh4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[nh4,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[nh4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[nh4,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[nh4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[nh4,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[nh4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[nh4,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/ammonium.jpg')
plt.close()

# Phosphate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phosphate [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[po4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[po4,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[po4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[po4,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[po4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[po4,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[po4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[po4,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[po4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[po4,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[po4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[po4,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[po4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[po4,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phosphate.jpg')
plt.close()

# Silicate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[sio4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[sio4,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[sio4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[sio4,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[sio4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[sio4,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[sio4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[sio4,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[sio4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[sio4,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[sio4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[sio4,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[sio4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[sio4,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/silicate.jpg')
plt.close()

# Bac1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Bac1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[bac1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[bac1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[bac1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[bac1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[bac1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[bac1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[bac1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Bac1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[bac1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[bac1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[bac1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[bac1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[bac1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[bac1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[bac1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Bac1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[bac1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[bac1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[bac1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[bac1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[bac1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[bac1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[bac1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1p.jpg')
plt.close()

# Phyto1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('phyto1 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto1l,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto1l,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto1l,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto1l,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto1l,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto1l,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto1l,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1l.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto1 Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1s,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto1s,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1s,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto1s,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1s,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto1s,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1s,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto1s,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1s,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto1s,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1s,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto1s,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1s,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto1s,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1s.jpg')
plt.close()

# Phyto2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto2c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto2c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto2c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto2c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto2c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto2c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto2c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto2n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto2n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto2n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto2n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto2n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto2n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto2n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto2p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto2p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto2p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto2p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto2p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto2p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto2p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto2 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto2l,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto2l,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto2l,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto2l,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto2l,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto2l,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto2l,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2l.jpg')
plt.close()

# Phyto3 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto3 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto3c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto3c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto3c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto3c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto3c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto3c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto3c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto3 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto3n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto3n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto3n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto3n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto3n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto3n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto3n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto3 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto3p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto3p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto3p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto3p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto3p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto3p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto3p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto3 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto3l,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto3l,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto3l,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto3l,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto3l,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto3l,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto3l,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3l.jpg')
plt.close()

# Phyto4 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto4 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto4c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto4c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto4c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto4c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto4c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto4c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto4c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto4 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto4n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto4n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto4n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto4n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto4n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto4n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto4n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto4 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto4p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto4p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto4p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto4p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto4p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto4p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto4p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phyto4 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[phyto4l,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[phyto4l,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[phyto4l,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[phyto4l,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[phyto4l,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[phyto4l,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[phyto4l,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4l.jpg')
plt.close()

# MesoZoo1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1p.jpg')
plt.close()

# MesoZoo2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz2c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz2c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz2c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz2c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz2c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz2c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz2c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz2n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz2n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz2n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz2n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz2n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz2n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz2n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MesoZoo2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[mesoz2p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[mesoz2p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[mesoz2p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[mesoz2p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[mesoz2p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[mesoz2p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[mesoz2p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2p.jpg')
plt.close()

# MicroZoo1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1p.jpg')
plt.close()

# MicroZoo2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz2c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz2c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz2c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz2c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz2c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz2c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz2c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz2n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz2n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz2n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz2n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz2n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz2n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz2n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('MicroZoo2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[microz2p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[microz2p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[microz2p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[microz2p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[microz2p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[microz2p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[microz2p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2p.jpg')
plt.close()

# Dom1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dom1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[dom1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[dom1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[dom1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[dom1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[dom1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[dom1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[dom1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dom1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[dom1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[dom1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[dom1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[dom1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[dom1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[dom1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[dom1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dom1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[dom1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[dom1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[dom1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[dom1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[dom1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[dom1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[dom1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1p.jpg')
plt.close()

# Dom2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dom2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[dom2c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[dom2c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[dom2c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[dom2c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[dom2c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[dom2c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[dom2c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom2c.jpg')
plt.close()

# Dom3 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dom3 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom3c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[dom3c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom3c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[dom3c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom3c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[dom3c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom3c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[dom3c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom3c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[dom3c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom3c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[dom3c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom3c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[dom3c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom3c.jpg')
plt.close()

# Pom1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Pom1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[pom1c,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[pom1c,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[pom1c,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[pom1c,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[pom1c,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[pom1c,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[pom1c,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Pom1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[pom1n,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[pom1n,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[pom1n,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[pom1n,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[pom1n,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[pom1n,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[pom1n,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Pom1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[pom1p,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[pom1p,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[pom1p,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[pom1p,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[pom1p,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[pom1p,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[pom1p,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Pom1 Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1s,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, pypom_daily[pom1s,0,:], '-.k', label='pyPOM')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1s,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, pypom_daily[pom1s,24,:], '-.k', label='pyPOM')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1s,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, pypom_daily[pom1s,49,:], '-.k', label='pyPOM')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1s,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, pypom_daily[pom1s,74,:], '-.k', label='pyPOM')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1s,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, pypom_daily[pom1s,99,:], '-.k', label='pyPOM')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1s,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, pypom_daily[pom1s,124,:], '-.k', label='pyPOM')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1s,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, pypom_daily[pom1s,149,:], '-.k', label='pyPOM')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1s.jpg')
plt.close()