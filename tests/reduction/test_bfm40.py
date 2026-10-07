from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os

# ----------------------------------------------------------------------------------------------------
# Extract model results
# ----------------------------------------------------------------------------------------------------
# Full BFM50
path = os.getcwd() + '/tests/bfm56/data/concentration_bfm56.npz'
globe_full = np.load(path, allow_pickle=True)

globe_daily = globe_full["daily"]
globe_monthly = globe_full["monthly"]

globe_daily = globe_daily[:,:,:1800]
globe_monthly = globe_monthly[:,:,:60]

# Tracer Indices
path = os.getcwd() + '/tests/bfm56/data/tracer_indices_bfm56.npz'
indices = np.load(path)
tracer_indices_full = {}
for file in indices.files:
    tracer_indices_full[file] = list(indices[file])

o2_full = globe_monthly[tracer_indices_full["o2"][0]]
no3_full = globe_monthly[tracer_indices_full["no3"][0]]
po4_full = globe_monthly[tracer_indices_full["po4"][0]]
chl_full = globe_monthly[tracer_indices_full["phyto1"][3]] + globe_monthly[tracer_indices_full["phyto2"][3]] + globe_monthly[tracer_indices_full["phyto3"][3]] + globe_monthly[tracer_indices_full["phyto4"][3]]
dic_full = globe_monthly[tracer_indices_full["co2"][0]]
pon_full = globe_monthly[tracer_indices_full["phyto1"][1]] + globe_monthly[tracer_indices_full["phyto2"][1]] + globe_monthly[tracer_indices_full["phyto3"][1]] + globe_monthly[tracer_indices_full["phyto4"][1]] \
            + globe_monthly[tracer_indices_full["mesozoo1"][1]] + globe_monthly[tracer_indices_full["mesozoo2"][1]] + globe_monthly[tracer_indices_full["microzoo1"][1]] + globe_monthly[tracer_indices_full["microzoo2"][1]] \
            + globe_monthly[tracer_indices_full["dom1"][1]] + globe_monthly[tracer_indices_full["pom1"][1]] 

conc_full = np.zeros((6,o2_full.shape[0],60))
conc_full[0] = chl_full[:,:60]
conc_full[1] = o2_full[:,:60]
conc_full[2] = no3_full[:,:60]
conc_full[3] = po4_full[:,:60]
conc_full[4] = pon_full[:,:60]
conc_full[5] = dic_full[:,:60]


# Reduced bfm40
path = os.getcwd() + '/tests/reduction/data/concentration_bfm40.npz'
bfm40 = np.load(path, allow_pickle=True)

bfm40_daily = bfm40["daily"]
bfm40_monthly = bfm40["monthly"]

path = os.getcwd() + '/tests/reduction/data/tracer_indices_bfm40.npz'
indices = np.load(path)
tracer_indices_bfm40 = {}
for file in indices.files:
    tracer_indices_bfm40[file] = list(indices[file])

o2_bfm40 = bfm40_monthly[tracer_indices_bfm40["o2"][0]]
no3_bfm40 = bfm40_monthly[tracer_indices_bfm40["no3"][0]]
po4_bfm40 = bfm40_monthly[tracer_indices_bfm40["po4"][0]]
chl_bfm40 = bfm40_monthly[tracer_indices_bfm40["phyto1"][3]] + bfm40_monthly[tracer_indices_bfm40["phyto2"][3]] + bfm40_monthly[tracer_indices_bfm40["phyto3"][3]] + bfm40_monthly[tracer_indices_bfm40["phyto4"][3]]
dic_bfm40 = bfm40_monthly[tracer_indices_bfm40["co2"][0]]
pon_bfm40 = bfm40_monthly[tracer_indices_bfm40["phyto1"][1]] + bfm40_monthly[tracer_indices_bfm40["phyto2"][1]] + bfm40_monthly[tracer_indices_bfm40["phyto3"][1]] + bfm40_monthly[tracer_indices_bfm40["phyto4"][1]] \
            + bfm40_monthly[tracer_indices_bfm40["microzoo1"][1]] + bfm40_monthly[tracer_indices_bfm40["microzoo2"][1]] + bfm40_monthly[tracer_indices_bfm40["dom1"][1]] + bfm40_monthly[tracer_indices_bfm40["pom1"][1]]

conc_bfm40 = np.zeros((6,o2_bfm40.shape[0],60))
conc_bfm40[0] = chl_bfm40[:,:60]
conc_bfm40[1] = o2_bfm40[:,:60]
conc_bfm40[2] = no3_bfm40[:,:60]
conc_bfm40[3] = po4_bfm40[:,:60]
conc_bfm40[4] = pon_bfm40[:,:60]
conc_bfm40[5] = dic_bfm40[:,:60]


# ----------------------------------------------------------------------------------------------------
# Plot model results
# ----------------------------------------------------------------------------------------------------
# title_globe = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phosphate','(e) PON','(f) NPP','(g) DIC']
# title_bfm56 = ['(h) Chl-a','(i) Oxygen','(j) Nitrate','(k) Phosphate','(l) PON','(m) NPP','(n) DIC']
title_full = ['(a)','(b)','(c)','(d)','(e)','(f)']
title_bfm40 = ['(g)','(h)','(i)','(j)','(k)','(l)']
# ---------------------------------------------------------------------------------------------------------------------------------
# Colorbar Limits
clow   = [0,180,0,0,0.1,30]
chigh  = [0.225,235,2.5,0.075,0.405,200]

fig,axes = plt.subplots(4,3,figsize=[16,15])
for i in range(0,6):
    plt.subplot(4,3,i+1)
    plt.imshow(conc_full[i,:,12:24],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    if i%3 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_full[i],fontsize=20)
    plt.clim(clow[i],chigh[i])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)

for i in range(6,12):
    plt.subplot(4,3,i+1)
    plt.imshow(conc_bfm40[i-6,:,12:24],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.yticks([0,50,100,150])
    if i%3 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_bfm40[i-6],fontsize=20)
    plt.clim(clow[i-6],chigh[i-6]) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)   

plt.tight_layout(h_pad=0.75, w_pad=0.75)

fig_name = os.getcwd() + '/tests/reduction/figures/bfm40_field_plots.jpg'
plt.savefig(fig_name)