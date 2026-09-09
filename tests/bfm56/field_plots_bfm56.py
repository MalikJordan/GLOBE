from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os

# Extract GLOBE data
globe_path = os.getcwd() + '/tests/bfm56/data/concentration_bfm56.npz'
model = np.load(globe_path, allow_pickle=True)
daily = model["daily"]
monthly = model["monthly"]

npp_path = os.getcwd() + '/tests/bfm56/data/npp_bfm56.npz'
npp = np.load(npp_path, allow_pickle=True)
npp_daily = npp["daily"]
npp_monthly = npp["monthly"]

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

# Extract fields of interest
globe_daily = np.zeros((7,daily.shape[1],daily.shape[2]))
globe_daily[0,:,:] = daily[phyto1l,:,:] + daily[phyto2l,:,:] + daily[phyto3l,:,:] + daily[phyto4l,:,:]  # Chlorophyll
globe_daily[1,:,:] = daily[o2,:,:]  # Oxygen
globe_daily[2,:,:] = daily[no3,:,:] # Nitrate
globe_daily[3,:,:] = daily[po4,:,:] # Phosphate
globe_daily[4,:,:] = daily[phyto1n,:,:] + daily[phyto2n,:,:] + daily[phyto3n,:,:] + daily[phyto4n,:,:] \
                   + daily[mesoz1n,:,:] + daily[mesoz2n,:,:] + daily[microz1n,:,:] + daily[microz2n,:,:] \
                   + daily[dom1n,:,:] + daily[pom1n,:,:]    # Particulate Organic Nitrogen
globe_daily[5,:,:] = npp_daily/12
globe_daily[6,:,:] = daily[co2,:,:]

globe_monthly = np.zeros((7,monthly.shape[1],monthly.shape[2]))
globe_monthly[0,:,:] = monthly[phyto1l,:,:] + monthly[phyto2l,:,:] + monthly[phyto3l,:,:] + monthly[phyto4l,:,:]  # Chlorophyll
globe_monthly[1,:,:] = monthly[o2,:,:]  # Oxygen
globe_monthly[2,:,:] = monthly[no3,:,:] # Nitrate
globe_monthly[3,:,:] = monthly[po4,:,:] # Phosphate
globe_monthly[4,:,:] = monthly[phyto1n,:,:] + monthly[phyto2n,:,:] + monthly[phyto3n,:,:] + monthly[phyto4n,:,:] \
                   + monthly[mesoz1n,:,:] + monthly[mesoz2n,:,:] + monthly[microz1n,:,:] + monthly[microz2n,:,:] \
                   + monthly[dom1n,:,:] + monthly[pom1n,:,:]    # Particulate Organic Nitrogen
globe_monthly[5,:,:] = npp_monthly/12
globe_monthly[6,:,:] = monthly[co2,:,:]


# Extract BFM data
bfm56_path = os.getcwd() + '/tests/bfm56/data/bfm56_pom1d.nc'
variables = nc.Dataset(bfm56_path)
variables = variables.variables

# Extract fields of interest
chlorophyll = (variables['Chla'][:])
oxygen = variables['O2o'][:]
nitrate = variables['N3n'][:]
phosphate = variables['N1p'][:]
pon = variables['P1n'][:] + variables['P2n'][:] + variables['P3n'][:] + variables['P4n'][:] + variables['Z3n'][:] + variables['Z4n'][:] + variables['Z5n'][:] + variables['Z6n'][:] \
    + variables['R1n'][:] + variables['R6n'][:]
# pon = variables['R6n'][:] + variables['P1n'][:] + variables['P2n'][:] + variables['P3n'][:] + variables['P4n'][:]
production = (variables['ruPTc'][:] - variables['resPP'][:] - variables['resZT'][:])/12
dic = (variables['DIC'][:])*(variables['ERHO'][:])*(12/1000)

# Write as array
chlorophyll = np.asarray(chlorophyll).transpose()
oxygen = np.asarray(oxygen).transpose()
nitrate = np.asarray(nitrate).transpose()
phosphate = np.asarray(phosphate).transpose()
pon = np.asarray(pon).transpose()
production = np.asarray(production).transpose()
dic = np.asarray(dic).transpose()

bfm56_daily = np.zeros((7,chlorophyll.shape[0],chlorophyll.shape[1]))
bfm56_daily[0,:,:] = chlorophyll
bfm56_daily[1,:,:] = oxygen
bfm56_daily[2,:,:] = nitrate
bfm56_daily[3,:,:] = phosphate
bfm56_daily[4,:,:] = pon
bfm56_daily[5,:,:] = production
bfm56_daily[6,:,:] = dic

# Monthly averages
bfm56_monthly = np.zeros((7,150,120))
for field in range(0,7):
    for year in range(0,10):
        for month in range(0,12):
            for day in range(0,30):
                bfm56_monthly[field,:,month + (year*12)] = bfm56_monthly[field,:,month + (year*12)] + bfm56_daily[field,:,(day + (month*30) + (year*360))]
bfm56_monthly = bfm56_monthly/30

bfm56_daily = bfm56_daily[:,:,:1800]
bfm56_monthly = bfm56_monthly[:,:,:60]

# Extract year 2 of data
globe_monthly = globe_monthly[:,:,12:24]
bfm56_monthly = bfm56_monthly[:,:,12:24]

# ---------------------------------------------------------------------------------------------------------------------------------
# Plot Style
plt.rc('font', family='serif', size=20)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=20, linewidth=1)
# ---------------------------------------------------------------------------------------------------------------------------------
# Legend Default
plt.rc('legend', framealpha=1.0, facecolor='white', frameon=True, edgecolor='black')
# ---------------------------------------------------------------------------------------------------------------------------------
# Plot Colors
bmap = brewer2mpl.get_map('Paired', 'qualitative', 10)
colors = bmap.mpl_colors
# ---------------------------------------------------------------------------------------------------------------------------------
# Titles
title_globe = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phosphate','(e) PON','(f) NPP','(g) DIC']
title_bfm56 = ['(h) Chl-a','(i) Oxygen','(j) Nitrate','(k) Phosphate','(l) PON','(m) NPP','(n) DIC']
title = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phosphate','(e) Chl-a','(f) Oxygen','(g) Nitrate','(h) Phoshate','(i) PON','(j) NPP','(k) DIC','(l) PON','(m) NPP','(n) DIC']
# ---------------------------------------------------------------------------------------------------------------------------------
# Colorbar Limits
clow   = [0,180,0,0,0.1,0,30]
chigh  = [0.225,235,2.5,0.075,0.405,2.0,170]
# ---------------------------------------------------------------------------------------------------------------------------------
# Field Plots  
    
fig,axes = plt.subplots(4,4,figsize=[16,15])
for i in range(0,7):
    plt.subplot(4,4,i+1)
    plt.imshow(globe_monthly[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
    plt.xlabel('Month',fontsize=14)
    if i%4 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_globe[i],fontsize=20)
    plt.clim(clow[i],chigh[i])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)

for i in range(7,14):
    plt.subplot(4,4,i+2)
    plt.imshow(bfm56_monthly[i-7,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
    plt.xlabel('Month',fontsize=14)
    plt.yticks([0,50,100,150])
    if i%4 == 3:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_bfm56[i-7],fontsize=20)
    plt.clim(clow[i-7],chigh[i-7]) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)   

fig.delaxes(axes[1,3])
fig.delaxes(axes[3,3])

plt.tight_layout(h_pad=0.75, w_pad=0.75)

fig_name = os.getcwd() + '/tests/bfm56/figures/bfm56_field_plots.jpg'
plt.savefig(fig_name)