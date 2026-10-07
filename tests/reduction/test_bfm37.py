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


# Reduced BFM37
path = os.getcwd() + '/tests/reduction/data/concentration_bfm37.npz'
bfm37 = np.load(path, allow_pickle=True)

bfm37_daily = bfm37["daily"]
bfm37_monthly = bfm37["monthly"]

path = os.getcwd() + '/tests/reduction/data/tracer_indices_bfm37.npz'
indices = np.load(path)
tracer_indices_bfm37 = {}
for file in indices.files:
    tracer_indices_bfm37[file] = list(indices[file])

o2_bfm37 = bfm37_monthly[tracer_indices_bfm37["o2"][0]]
no3_bfm37 = bfm37_monthly[tracer_indices_bfm37["no3"][0]]
po4_bfm37 = bfm37_monthly[tracer_indices_bfm37["po4"][0]]
chl_bfm37 = bfm37_monthly[tracer_indices_bfm37["phyto1"][3]] + bfm37_monthly[tracer_indices_bfm37["phyto2"][3]] + bfm37_monthly[tracer_indices_bfm37["phyto3"][3]] + bfm37_monthly[tracer_indices_bfm37["phyto4"][3]]
# dic_bfm37 = bfm37_monthly[tracer_indices_bfm37["co2"][0]]
pon_bfm37 = bfm37_monthly[tracer_indices_bfm37["phyto1"][1]] + bfm37_monthly[tracer_indices_bfm37["phyto2"][1]] + bfm37_monthly[tracer_indices_bfm37["phyto3"][1]] + bfm37_monthly[tracer_indices_bfm37["phyto4"][1]] \
            + bfm37_monthly[tracer_indices_bfm37["microzoo1"][1]] + bfm37_monthly[tracer_indices_bfm37["microzoo2"][1]] + bfm37_monthly[tracer_indices_bfm37["dom1"][1]] + bfm37_monthly[tracer_indices_bfm37["pom1"][1]]

conc_bfm37 = np.zeros((6,o2_bfm37.shape[0],60))
conc_bfm37[0] = chl_bfm37[:,:60]
conc_bfm37[1] = o2_bfm37[:,:60]
conc_bfm37[2] = no3_bfm37[:,:60]
conc_bfm37[3] = po4_bfm37[:,:60]
conc_bfm37[4] = pon_bfm37[:,:60]
# conc_bfm37[5] = dic_bfm37[:,:60]

# Reduced BFM37

# bfm37_daily = bfm37["daily"]

# path = os.getcwd() + '/tests/reduction/data/tracer_indices_bfm37.npz'
# indices = np.load(path)
# tracer_indices_bfm37 = {}
# for file in indices.files:
#     tracer_indices_bfm37[file] = list(indices[file])


# chl_full = globe_daily[13] + globe_daily[18] + globe_daily[22] + globe_daily[26]
# chl_bfm37 = bfm37_daily[10] + bfm37_daily[15] + bfm37_daily[19] + bfm37_daily[23]



# path_full = os.getcwd() + '/bfm50_ddt_iter0.npz'
# path_bfm37 = os.getcwd() + '/bfm37_ddt_iter0.npz'

# d_dt50 = np.load(path_full,allow_pickle=True)["d_dt"]
# d_dt37 = np.load(path_bfm37,allow_pickle=True)["d_dt"]

# tracers37 = [0,1,2,3,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,33,34,35,36,37,38,39,40,41,44,45,46,47]
# d_dt50 = d_dt50[tracers37,:]
# dif = d_dt50 - d_dt37
# x=1



# ----------------------------------------------------------------------------------------------------
# Plot model results
# ----------------------------------------------------------------------------------------------------
# title_globe = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phosphate','(e) PON','(f) NPP','(g) DIC']
# title_bfm56 = ['(h) Chl-a','(i) Oxygen','(j) Nitrate','(k) Phosphate','(l) PON','(m) NPP','(n) DIC']
title_full = ['(a)','(b)','(c)','(d)','(e)','(f)']
title_bfm37 = ['(g)','(h)','(i)','(j)','(k)','(l)']
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
    plt.imshow(conc_bfm37[i-6,:,12:24],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.yticks([0,50,100,150])
    if i%3 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_bfm37[i-6],fontsize=20)
    plt.clim(clow[i-6],chigh[i-6]) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)   

plt.tight_layout(h_pad=0.75, w_pad=0.75)

fig_name = os.getcwd() + '/tests/reduction/figures/bfm37_field_plots.jpg'
plt.savefig(fig_name)