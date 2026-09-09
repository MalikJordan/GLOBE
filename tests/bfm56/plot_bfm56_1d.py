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
path = os.getcwd() + '/tests/bfm56/data/concentration_bfm56.npz'
globe_bfm = np.load(path, allow_pickle=True)

globe_daily = globe_bfm["daily"]
globe_monthly = globe_bfm["monthly"]

globe_daily = globe_daily[:,:,:1800]
globe_monthly = globe_monthly[:,:,:60]

# Tracer Indices
path = os.getcwd() + '/tests/bfm56/data/tracer_indices_bfm56.npz'
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
path = os.getcwd() + '/tests/bfm56/data/bfm56_pom1d.nc'
variables = nc.Dataset(path)
variables = variables.variables

o2o = np.asarray(variables['O2o'][:]).transpose()   # oxygen
n1p = np.asarray(variables['N1p'][:]).transpose()   # phosphate
n3n = np.asarray(variables['N3n'][:]).transpose()   # nitrate
n4n = np.asarray(variables['N4n'][:]).transpose()   # ammonium
o4n = np.asarray(variables['O4n'][:]).transpose()   # nitrate sink
n5s  = np.asarray(variables['N5s'][:]).transpose()  # silicate
n6r  = np.asarray(variables['N6r'][:]).transpose()  # reduction equivalents
b1c = np.asarray(variables['B1c'][:]).transpose()   # bacteria (anerobic and aerobic)
b1n = np.asarray(variables['B1n'][:]).transpose()
b1p = np.asarray(variables['B1p'][:]).transpose()
p1c = np.asarray(variables['P1c'][:]).transpose()   # diatoms
p1n = np.asarray(variables['P1n'][:]).transpose()
p1p = np.asarray(variables['P1p'][:]).transpose()
p1l = np.asarray(variables['P1l'][:]).transpose()
p1s = np.asarray(variables['P1s'][:]).transpose()
p2c = np.asarray(variables['P2c'][:]).transpose()   # flagellates
p2n = np.asarray(variables['P2n'][:]).transpose()
p2p = np.asarray(variables['P2p'][:]).transpose()
p2l = np.asarray(variables['P2l'][:]).transpose()
p3c = np.asarray(variables['P3c'][:]).transpose()   # picophytoplankton
p3n = np.asarray(variables['P3n'][:]).transpose()
p3p = np.asarray(variables['P3p'][:]).transpose()
p3l = np.asarray(variables['P3l'][:]).transpose()
p4c = np.asarray(variables['P4c'][:]).transpose()   # large phytoplankton
p4n = np.asarray(variables['P4n'][:]).transpose()
p4p = np.asarray(variables['P4p'][:]).transpose()
p4l = np.asarray(variables['P4l'][:]).transpose()
z3c = np.asarray(variables['Z3c'][:]).transpose()   # carnivorous mesozooplankton
z3n = np.asarray(variables['Z3n'][:]).transpose()
z3p = np.asarray(variables['Z3p'][:]).transpose()
z4c = np.asarray(variables['Z4c'][:]).transpose()   # omnivorous mesozooplankton
z4n = np.asarray(variables['Z4n'][:]).transpose()
z4p = np.asarray(variables['Z4p'][:]).transpose()
z5c = np.asarray(variables['Z5c'][:]).transpose()   # microzooplankton
z5n = np.asarray(variables['Z5n'][:]).transpose()
z5p = np.asarray(variables['Z5p'][:]).transpose()
z6c = np.asarray(variables['Z6c'][:]).transpose()   # heterotrophic nanoflagellates
z6n = np.asarray(variables['Z6n'][:]).transpose()
z6p = np.asarray(variables['Z6p'][:]).transpose()
r1c = np.asarray(variables['R1c'][:]).transpose()   # labile dissolved organic matter
r1n = np.asarray(variables['R1n'][:]).transpose()
r1p = np.asarray(variables['R1p'][:]).transpose()
r2c = np.asarray(variables['R2c'][:]).transpose()   # semi-labile dissolved organic carbon
r3c = np.asarray(variables['R3c'][:]).transpose()   # semi-refractory dissolved organic carbon
r6c = np.asarray(variables['R6c'][:]).transpose()   # particulate organic matter
r6n = np.asarray(variables['R6n'][:]).transpose()
r6p = np.asarray(variables['R6p'][:]).transpose()
r6s = np.asarray(variables['R6s'][:]).transpose()
o3c = np.asarray(variables['O3c'][:]).transpose()   # dissolved inorganic carbon
o3h = np.asarray(variables['O3h'][:]).transpose()   # total alkalinity




bfm56_daily = np.zeros((50,o2o.shape[0],o2o.shape[1]))
bfm56_daily[0,:,:] = o2o
bfm56_daily[1,:,:] = n1p
bfm56_daily[2,:,:] = n3n
bfm56_daily[3,:,:] = n4n
bfm56_daily[4,:,:] = o4n
bfm56_daily[5,:,:] = n5s
bfm56_daily[6,:,:] = n6r
bfm56_daily[7,:,:] = b1c
bfm56_daily[8,:,:] = b1n
bfm56_daily[9,:,:] = b1p
bfm56_daily[10,:,:] = p1c
bfm56_daily[11,:,:] = p1n
bfm56_daily[12,:,:] = p1p
bfm56_daily[13,:,:] = p1l
bfm56_daily[14,:,:] = p1s
bfm56_daily[15,:,:] = p2c
bfm56_daily[16,:,:] = p2n
bfm56_daily[17,:,:] = p2p
bfm56_daily[18,:,:] = p2l
bfm56_daily[19,:,:] = p3c
bfm56_daily[20,:,:] = p3n
bfm56_daily[21,:,:] = p3p
bfm56_daily[22,:,:] = p3l
bfm56_daily[23,:,:] = p4c
bfm56_daily[24,:,:] = p4n
bfm56_daily[25,:,:] = p4p
bfm56_daily[26,:,:] = p4l
bfm56_daily[27,:,:] = z3c
bfm56_daily[28,:,:] = z3n
bfm56_daily[29,:,:] = z3p
bfm56_daily[30,:,:] = z4c
bfm56_daily[31,:,:] = z4n
bfm56_daily[32,:,:] = z4p
bfm56_daily[33,:,:] = z5c
bfm56_daily[34,:,:] = z5n
bfm56_daily[35,:,:] = z5p
bfm56_daily[36,:,:] = z6c
bfm56_daily[37,:,:] = z6n
bfm56_daily[38,:,:] = z6p
bfm56_daily[39,:,:] = r1c
bfm56_daily[40,:,:] = r1n
bfm56_daily[41,:,:] = r1p
bfm56_daily[42,:,:] = r2c
bfm56_daily[43,:,:] = r3c
bfm56_daily[44,:,:] = r6c
bfm56_daily[45,:,:] = r6n
bfm56_daily[46,:,:] = r6p
bfm56_daily[47,:,:] = r6s
bfm56_daily[48,:,:] = o3c
bfm56_daily[49,:,:] = o3h

bfm56_monthly = np.zeros((50,150,120))
for spec in range(0,50):
    for year in range(0,10):
        for month in range(0,12):
            for day in range(0,30):
                bfm56_monthly[spec,:,month + (year*12)] = bfm56_monthly[spec,:,month + (year*12)] + bfm56_daily[spec,:,(day + (month*30) + (year*360))]
bfm56_monthly = bfm56_monthly/30

bfm56_daily = bfm56_daily[:,:,:1800]
bfm56_monthly = bfm56_monthly[:,:,:60]

# ----------------------------------------------------------------------------------------------------
# Field Plots
# ----------------------------------------------------------------------------------------------------
pon_globe = globe_monthly[11] + globe_monthly[16] + globe_monthly[20] + globe_monthly[24] + globe_monthly[28] + globe_monthly[31] + globe_monthly[34] + globe_monthly[37] + globe_monthly[40] + globe_monthly[45]
pon_bfm56 = bfm56_monthly[11] + bfm56_monthly[16] + bfm56_monthly[20] + bfm56_monthly[24] + bfm56_monthly[28] + bfm56_monthly[31] + bfm56_monthly[34] + bfm56_monthly[37] + bfm56_monthly[40] + bfm56_monthly[45]

chl_globe = globe_monthly[13] + globe_monthly[18] + globe_monthly[22] + globe_monthly[26]
chl_bfm56 = globe_monthly[13] + globe_monthly[18] + globe_monthly[22] + globe_monthly[26]

# Colorbar Limits
clow   = [0,180,0,0,0.1,30]
chigh  = [0.225,235,2.5,0.075,0.405,170]

title_globe = ['GLOBE - Chl-a','GLOBE - Oxygen','GLOBE - Nitrate','GLOBE - Phosphate','GLOBE - PON','GLOBE - DIC']
title_bfm56 = ['BFM56 - Chl-a','BFM56 - Oxygen','BFM56 - Nitrate','BFM56 - Phosphate','BFM56 - PON','BFM56 - DIC']

# fig,ax = plt.subplots(2,6,figsize=[20,10])
# for i in range(7):

#     ax[0,i].imshow(chl_globe[:,12:24],extent=[0,12,150,0],aspect='auto')

# ----------------------------------------------------------------------------------------------------
# Line Plots
# ----------------------------------------------------------------------------------------------------
days = np.linspace(0,1799,1800)
xticks = [0,360,720,1080,1440]
xlabels_top = ['','','','','']
xlabels_bottom = ['1','2','3','4','5']

# Oxygen -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Oxygen [mmol $O_{2}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[o2,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[o2,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[o2,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[o2,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[o2,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[o2,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[o2,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[o2,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[o2,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[o2,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[o2,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[o2,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[o2,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[o2,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/oxygen.jpg')
plt.close()

# Dissolved Inorganic Carbon -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dissolved Inorganic Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[co2,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[co2,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[co2,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[co2,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[co2,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[co2,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[co2,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[co2,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[co2,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[co2,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[co2,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[co2,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[co2,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[co2,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/co2.jpg')
plt.close()

# Nitrate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Nitrate [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[no3,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[no3,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[no3,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[no3,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[no3,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[no3,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[no3,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[no3,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[no3,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[no3,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[no3,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[no3,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[no3,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[no3,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/nitrate.jpg')
plt.close()

# Ammonium -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Ammonium [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[nh4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[nh4,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[nh4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[nh4,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[nh4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[nh4,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[nh4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[nh4,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[nh4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[nh4,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[nh4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[nh4,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[nh4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[nh4,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/ammonium.jpg')
plt.close()

# Phosphate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phosphate [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[po4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[po4,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[po4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[po4,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[po4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[po4,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[po4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[po4,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[po4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[po4,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[po4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[po4,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[po4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[po4,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phosphate.jpg')
plt.close()

# Silicate -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[sio4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[sio4,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[sio4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[sio4,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[sio4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[sio4,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[sio4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[sio4,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[sio4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[sio4,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[sio4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[sio4,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[sio4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[sio4,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/silicate.jpg')
plt.close()

# Bac1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Bac1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[bac1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[bac1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[bac1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[bac1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[bac1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[bac1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[bac1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Bac1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[bac1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[bac1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[bac1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[bac1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[bac1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[bac1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[bac1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Bac1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[bac1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[bac1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[bac1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[bac1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[bac1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[bac1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[bac1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[bac1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[bac1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[bac1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[bac1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[bac1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[bac1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[bac1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/bac1p.jpg')
plt.close()

# Phyto1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('phyto1 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto1l,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto1l,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto1l,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto1l,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto1l,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto1l,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto1l,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1l.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto1 Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto1s,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto1s,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto1s,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto1s,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto1s,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto1s,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto1s,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto1s,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto1s,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto1s,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto1s,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto1s,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto1s,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto1s,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto1s.jpg')
plt.close()

# Phyto2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto2c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto2c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto2c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto2c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto2c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto2c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto2c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto2n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto2n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto2n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto2n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto2n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto2n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto2n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto2p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto2p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto2p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto2p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto2p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto2p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto2p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto2 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto2l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto2l,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto2l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto2l,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto2l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto2l,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto2l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto2l,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto2l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto2l,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto2l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto2l,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto2l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto2l,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto2l.jpg')
plt.close()

# Phyto3 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto3 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto3c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto3c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto3c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto3c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto3c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto3c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto3c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto3 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto3n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto3n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto3n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto3n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto3n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto3n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto3n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto3 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto3p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto3p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto3p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto3p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto3p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto3p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto3p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto3 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto3l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto3l,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto3l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto3l,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto3l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto3l,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto3l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto3l,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto3l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto3l,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto3l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto3l,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto3l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto3l,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto3l.jpg')
plt.close()

# Phyto4 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto4 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto4c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto4c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto4c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto4c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto4c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto4c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto4c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto4 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto4n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto4n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto4n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto4n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto4n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto4n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto4n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto4 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto4p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto4p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto4p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto4p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto4p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto4p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto4p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Phyto4 Chlorophyll-a [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[phyto4l,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[phyto4l,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[phyto4l,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[phyto4l,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[phyto4l,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[phyto4l,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[phyto4l,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[phyto4l,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[phyto4l,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[phyto4l,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[phyto4l,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[phyto4l,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[phyto4l,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[phyto4l,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/phyto4l.jpg')
plt.close()

# MesoZoo1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz1p.jpg')
plt.close()

# MesoZoo2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz2c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz2c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz2c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz2c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz2c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz2c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz2c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz2n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz2n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz2n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz2n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz2n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz2n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz2n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MesoZoo2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[mesoz2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[mesoz2p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[mesoz2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[mesoz2p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[mesoz2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[mesoz2p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[mesoz2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[mesoz2p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[mesoz2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[mesoz2p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[mesoz2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[mesoz2p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[mesoz2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[mesoz2p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/mesoz2p.jpg')
plt.close()

# MicroZoo1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz1p.jpg')
plt.close()

# MicroZoo2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz2c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz2c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz2c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz2c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz2c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz2c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz2c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo2 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz2n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz2n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz2n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz2n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz2n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz2n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz2n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('MicroZoo2 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[microz2p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[microz2p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[microz2p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[microz2p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[microz2p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[microz2p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[microz2p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[microz2p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[microz2p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[microz2p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[microz2p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[microz2p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[microz2p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[microz2p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/microz2p.jpg')
plt.close()

# Dom1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dom1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[dom1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[dom1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[dom1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[dom1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[dom1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[dom1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[dom1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dom1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[dom1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[dom1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[dom1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[dom1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[dom1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[dom1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[dom1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dom1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[dom1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[dom1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[dom1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[dom1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[dom1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[dom1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[dom1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom1p.jpg')
plt.close()

# Dom2 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dom2 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom2c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[dom2c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom2c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[dom2c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom2c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[dom2c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom2c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[dom2c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom2c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[dom2c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom2c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[dom2c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom2c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[dom2c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom2c.jpg')
plt.close()

# Dom3 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Dom3 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[dom3c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[dom3c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[dom3c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[dom3c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[dom3c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[dom3c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[dom3c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[dom3c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[dom3c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[dom3c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[dom3c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[dom3c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[dom3c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[dom3c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/dom3c.jpg')
plt.close()

# Pom1 -----------------------------------------------------------------------------------------------------------
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Pom1 Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1c,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[pom1c,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1c,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[pom1c,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1c,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[pom1c,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1c,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[pom1c,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1c,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[pom1c,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1c,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[pom1c,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1c,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[pom1c,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1c.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Pom1 Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1n,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[pom1n,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1n,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[pom1n,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1n,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[pom1n,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1n,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[pom1n,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1n,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[pom1n,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1n,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[pom1n,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1n,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[pom1n,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1n.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Pom1 Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1p,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[pom1p,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1p,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[pom1p,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1p,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[pom1p,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1p,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[pom1p,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1p,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[pom1p,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1p,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[pom1p,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1p,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[pom1p,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1p.jpg')
plt.close()

fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
plt.suptitle('Pom1 Silicate [mmol $Si$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, globe_daily[pom1s,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, bfm56_daily[pom1s,0,:], '-.k', label='BFM56')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, globe_daily[pom1s,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, bfm56_daily[pom1s,24,:], '-.k', label='BFM56')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, globe_daily[pom1s,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, bfm56_daily[pom1s,49,:], '-.k', label='BFM56')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, globe_daily[pom1s,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, bfm56_daily[pom1s,74,:], '-.k', label='BFM56')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, globe_daily[pom1s,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, bfm56_daily[pom1s,99,:], '-.k', label='BFM56')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, globe_daily[pom1s,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, bfm56_daily[pom1s,124,:], '-.k', label='BFM56')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, globe_daily[pom1s,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, bfm56_daily[pom1s,149,:], '-.k', label='BFM56')

for i in range(4):
    ax[0,i].set_xlim(0,1800)
    ax[0,i].set_xticks(xticks)

    if i != 3:
        ax[1,i].set_xlim(0,1800)
        ax[1,i].set_xticks(xticks)
        ax[1,i].set_xticklabels(['1','2','3','4','5'])
        ax[1,i].set_xlabel('Year')

handles,labels = ax[0,0].get_legend_handles_labels()
ax[1,3].axis('off')
ax[1,3].legend(handles,labels,loc='center',frameon=True)

plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm56/figures/pom1s.jpg')
plt.close()