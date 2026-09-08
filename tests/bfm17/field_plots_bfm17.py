from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os

# Extract GLOBE data
globe_path = os.getcwd() + '/concentration_bfm17-5yr-0907.npz'
model = np.load(globe_path, allow_pickle=True)
daily = model["daily"]
monthly = model["monthly"]

npp_path = os.getcwd() + '/npp_bfm17-5yr-0907.npz'
npp = np.load(npp_path, allow_pickle=True)
npp_daily = npp["daily"]
npp_monthly = npp["monthly"]

# Tracer Indices
path = os.getcwd() + '/tracer_indices_bfm17-5yr.npz'
indices = np.load(path)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

o2 = tracer_indices["o2"][0]
no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
po4 = tracer_indices["po4"][0]

phyto1c = tracer_indices["phyto1"][0]
phyto1n = tracer_indices["phyto1"][1]
phyto1p = tracer_indices["phyto1"][2]
phyto1l = tracer_indices["phyto1"][3]

microz1c = tracer_indices["zoo1"][0]
microz1n = tracer_indices["zoo1"][1]
microz1p = tracer_indices["zoo1"][2]

dom1c = tracer_indices["dom1"][0]
dom1n = tracer_indices["dom1"][1]
dom1p = tracer_indices["dom1"][2]

pom1c = tracer_indices["pom1"][0]
pom1n = tracer_indices["pom1"][1]
pom1p = tracer_indices["pom1"][2]

# Extract fields of interest
globe_daily = np.zeros((7,daily.shape[1],daily.shape[2]))
globe_daily[0,:,:] = daily[phyto1l,:,:] # Chlorophyll
globe_daily[1,:,:] = daily[o2,:,:]  # Oxygen
globe_daily[2,:,:] = daily[no3,:,:] # Nitrate
globe_daily[3,:,:] = daily[po4,:,:] # Phosphate
globe_daily[4,:,:] = daily[phyto1n,:,:] + daily[microz1n,:,:] + daily[dom1n,:,:] + daily[pom1n,:,:]    # Particulate Organic Nitrogen
globe_daily[5,:,:] = npp_daily/12

globe_monthly = np.zeros((7,monthly.shape[1],monthly.shape[2]))
globe_monthly[0,:,:] = monthly[phyto1l,:,:] # Chlorophyll
globe_monthly[1,:,:] = monthly[o2,:,:]  # Oxygen
globe_monthly[2,:,:] = monthly[no3,:,:] # Nitrate
globe_monthly[3,:,:] = monthly[po4,:,:] # Phosphate
globe_monthly[4,:,:] = monthly[phyto1n,:,:] + monthly[microz1n,:,:] + monthly[dom1n,:,:] + monthly[pom1n,:,:]    # Particulate Organic Nitrogen
globe_monthly[5,:,:] = npp_monthly/12


# Extract BFM data
bfm17_path = os.getcwd() + '/tests/bfm17/data/bfm17_pom1d.nc'
variables = nc.Dataset(bfm17_path)
variables = variables.variables

# Extract fields of interest
chlorophyll = (variables['Chla'][:])
oxygen = variables['O2o'][:]
nitrate = variables['N3n'][:]
phosphate = variables['N1p'][:]
pon = variables['P2n'][:] + variables['Z5n'][:] + variables['R1n'][:] + variables['R6n'][:]
# pon = variables['R6n'][:] + variables['P1n'][:] + variables['P2n'][:] + variables['P3n'][:] + variables['P4n'][:]
production = (variables['ruPTc'][:] - variables['resPP'][:])/12
# production = (variables['ruPTc'][:] - variables['resPP'][:] - variables['resZT'][:])/12

# Write as array
chlorophyll = np.asarray(chlorophyll).transpose()
oxygen = np.asarray(oxygen).transpose()
nitrate = np.asarray(nitrate).transpose()
phosphate = np.asarray(phosphate).transpose()
pon = np.asarray(pon).transpose()
production = np.asarray(production).transpose()

bfm17_daily = np.zeros((6,chlorophyll.shape[0],chlorophyll.shape[1]))
bfm17_daily[0,:,:] = chlorophyll
bfm17_daily[1,:,:] = oxygen
bfm17_daily[2,:,:] = nitrate
bfm17_daily[3,:,:] = phosphate
bfm17_daily[4,:,:] = pon
bfm17_daily[5,:,:] = production

# Monthly averages
bfm17_monthly = np.zeros((6,150,120))
for field in range(0,6):
    for year in range(0,10):
        for month in range(0,12):
            for day in range(0,30):
                bfm17_monthly[field,:,month + (year*12)] = bfm17_monthly[field,:,month + (year*12)] + bfm17_daily[field,:,(day + (month*30) + (year*360))]
bfm17_monthly = bfm17_monthly/30

bfm17_daily = bfm17_daily[:,:,:1800]
bfm17_monthly = bfm17_monthly[:,:,:60]

# Extract year 2 of data
globe_monthly = globe_monthly[:,:,12:24]
bfm17_monthly = bfm17_monthly[:,:,12:24]

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
title_globe = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phosphate','(e) PON','(f) NPP']
title_bfm17 = ['(h) Chl-a','(i) Oxygen','(j) Nitrate','(k) Phosphate','(l) PON','(m) NPP']
# ---------------------------------------------------------------------------------------------------------------------------------
# Colorbar Limits
clow   = [0,180,0,0,0.1,0]
chigh  = [0.225,235,2.5,0.075,0.405,2.0]
# ---------------------------------------------------------------------------------------------------------------------------------
# Field Plots  
    
fig,axes = plt.subplots(4,3,figsize=[16,15])
for i in range(0,6):
    plt.subplot(4,3,i+1)
    plt.imshow(globe_monthly[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
    plt.xlabel('Month',fontsize=14)
    if i%3 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_globe[i],fontsize=20)
    plt.clim(clow[i],chigh[i])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)

for i in range(6,12):
    plt.subplot(4,3,i+1)
    plt.imshow(bfm17_monthly[i-6,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
    ax = plt.gca()
    plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
    plt.xlabel('Month',fontsize=14)
    plt.yticks([0,50,100,150])
    if i%3 == 0:
        plt.yticks([0,50,100,150])
        plt.ylabel('Depth (m)',fontsize=14)
    else:
        plt.yticks([0,50,100,150],[])
    plt.title(title_bfm17[i-6],fontsize=20)
    plt.clim(clow[i-6],chigh[i-6]) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)   

plt.tight_layout(h_pad=0.75, w_pad=0.75)

fig_name = os.getcwd() + '/tests/bfm17/figures/bfm17_field_plots.jpg'
plt.savefig(fig_name)