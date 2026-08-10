from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os
import colormaps as cmaps


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
# path = os.getcwd() + '/tests/bfm17/data/bfm17_pom1d_case44.nc'
path = os.getcwd() + '/tests/bfm17/data/bfm17_pom1d_dp64.nc'
variables = nc.Dataset(path)
variables = variables.variables

# chlorophyll = np.asarray(variables['Chla'][:]).transpose()
oxygen = np.asarray(variables['O2o'][:]).transpose()
nitrate = np.asarray(variables['N3n'][:]).transpose()
ammonium = np.asarray(variables['N4n'][:]).transpose()
phosphate = np.asarray(variables['N1p'][:]).transpose()
phytoc = np.asarray(variables['P2c'][:]).transpose()
phyton = np.asarray(variables['P2n'][:]).transpose()
phytop = np.asarray(variables['P2p'][:]).transpose()
phytol = np.asarray(variables['P2l'][:]).transpose()
zooc = np.asarray(variables['Z5c'][:]).transpose()
zoon = np.asarray(variables['Z5n'][:]).transpose()
zoop = np.asarray(variables['Z5p'][:]).transpose()
domc = np.asarray(variables['R1c'][:]).transpose()
domn = np.asarray(variables['R1n'][:]).transpose()
domp = np.asarray(variables['R1p'][:]).transpose()
pomc = np.asarray(variables['R6c'][:]).transpose()
pomn = np.asarray(variables['R6n'][:]).transpose()
pomp = np.asarray(variables['R6p'][:]).transpose()

data_fortran = np.zeros((17,oxygen.shape[0],oxygen.shape[1]))
data_fortran[0,:,:] = oxygen
data_fortran[1,:,:] = nitrate
data_fortran[2,:,:] = ammonium
data_fortran[3,:,:] = phosphate
data_fortran[4,:,:] = phytoc
data_fortran[5,:,:] = phyton
data_fortran[6,:,:] = phytop
data_fortran[7,:,:] = phytol
data_fortran[8,:,:] = zooc
data_fortran[9,:,:] = zoon
data_fortran[10,:,:] = zoop
data_fortran[11,:,:] = domc
data_fortran[12,:,:] = domn
data_fortran[13,:,:] = domp
data_fortran[14,:,:] = pomc
data_fortran[15,:,:] = pomn
data_fortran[16,:,:] = pomp

# avg_data_fortran = np.zeros((17,150,12))
# for spec in range(0,17):
#     for year in range(1,2):
#     # for year in range(14,15):
#         for month in range(0,12):
#             for day in range(0,30):
#                 avg_data_fortran[spec,:,month] = avg_data_fortran[spec,:,month] + data_fortran[spec,:,(day + (month*30) + (year*360))]
# avg_data_fortran = avg_data_fortran/30

# data_fortran = data_fortran[:,:,:1800]
# avg_data_fortran = avg_data_fortran[:,:,:60]

avg_data_fortran = np.zeros((17,150,12))
for spec in range(0,17):
    for year in range(0,1):
    # for year in range(14,15):
        for month in range(0,12):
            for day in range(0,30):
                avg_data_fortran[spec,:,month] = avg_data_fortran[spec,:,month] + data_fortran[spec,:,(day + (month*30) + (year*360))]
avg_data_fortran = avg_data_fortran/30

data_fortran = data_fortran[:,:,:360]
avg_data_fortran = avg_data_fortran[:,:,:12]


# for spec in range(0,17):
#     for year in range(1,2):
#     # for year in range(14,15):
#         for month in range(0,12):
#             for day in range(0,30):
#                 avg_data_fortran[spec,:,month+12] = avg_data_fortran[spec,:,month+12] + data_fortran[spec,:,(day + (month*30) + (year*360))]
# avg_data_fortran = avg_data_fortran/30


# path = os.getcwd() + '/tests/bfm17/data/concentration_bfm17-1d-old.npz'
path = os.getcwd() + '/concentration_bfm17-1d.npz'
# path = os.getcwd() + '/concentration_bfm56-5yr.npz'
model = np.load(path, allow_pickle=True)

bfm17_daily = model["daily"]
bfm17_monthly = model["monthly"]

# days = np.linspace(0,3599,3600)
# xticks = [0,360,720,1080,1440,1800,2160,2520,2880,3240]
# xlabel = ['1','2','3','4','5','6','7','8','9','10']

# bfm17_daily = bfm17_daily[:,:,:1800]
# bfm17_monthly = bfm17_monthly[:,:,:60]
# days = np.linspace(0,1799,1800)


days = np.linspace(0,359,360)
xticks = [15,45,75,105,135,165,195,225,255,285,315,345]
xlabel = ['J','','','A','','','J','','','O','','']

bfm17_daily = bfm17_daily[:,:,:360]
bfm17_monthly = bfm17_monthly[:,:,:12]
days = np.linspace(0,359,360)

# ----------------------------------------------------------------------------------------------------
# Line Plots
# ----------------------------------------------------------------------------------------------------
# Oxygen
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Oxygen [mmol $O_{2}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[0,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[0,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[0,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[0,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[0,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[0,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[0,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[0,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[0,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[0,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[0,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[0,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[0,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[0,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/oxygen.jpg')

# Nitrate
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Nitrate [mmol $NO_{3}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[1,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[1,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[1,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[1,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[1,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[1,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[1,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[1,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[1,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[1,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[1,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[1,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[1,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[1,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/nitrate.jpg')

# Ammonium
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Ammonium [mmol $NH_{4}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[2,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[2,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[2,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[2,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[2,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[2,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[2,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[2,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[2,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[2,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[2,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[2,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[2,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[2,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/ammonium.jpg')

# Phosphate
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phosphate [mmol $PO_{4}$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[3,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[3,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[3,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[3,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[3,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[3,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[3,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[3,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[3,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[3,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[3,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[3,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[3,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[3,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phosphate.jpg')

# Phyto Carbon
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phytoplankton Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[4,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[4,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[4,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[4,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[4,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[4,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[4,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[4,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[4,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[4,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[4,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[4,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[4,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[4,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phytoc.jpg')

# Phyto Nitrogen
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phytoplankton Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[5,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[5,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[5,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[5,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[5,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[5,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[5,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[5,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[5,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[5,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[5,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[5,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[5,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[5,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phyton.jpg')

# Phyto Phosphorus
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phytoplankton Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[6,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[6,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[6,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[6,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[6,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[6,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[6,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[6,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[6,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[6,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[6,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[6,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[6,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[6,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phytop.jpg')

# Phyto Chlorophyll
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Phytoplankton Chlorophyll [mg $Chl-a$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[7,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[7,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[7,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[7,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[7,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[7,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[7,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[7,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[7,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[7,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[7,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[7,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[7,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[7,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phytol.jpg')

# Zoo Carbon
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Zooplankton Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[8,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[8,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[8,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[8,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[8,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[8,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[8,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[8,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[8,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[8,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[8,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[8,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[8,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[8,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/zooc.jpg')

# Zoo Nitrogen
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Zooplankton Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[9,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[9,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[9,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[9,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[9,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[9,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[9,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[9,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[9,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[9,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[9,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[9,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[9,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[9,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/zoon.jpg')

# Zoo Phosphorus
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Zooplankton Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[10,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[10,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[10,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[10,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[10,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[10,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[10,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[10,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[10,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[10,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[10,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[10,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[10,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[10,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/zoop.jpg')

# DOM Carbon
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dissolved Organic Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[11,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[11,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[11,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[11,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[11,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[11,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[11,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[11,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[11,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[11,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[11,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[11,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[11,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[11,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/domc.jpg')

# DOM Nitrogen
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dissolved Organic Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[12,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[12,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[12,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[12,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[12,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[12,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[12,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[12,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[12,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[12,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[12,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[12,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[12,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[12,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/domn.jpg')

# DOM Phosphorus
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Dissolved Organic Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[13,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[13,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[13,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[13,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[13,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[13,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[13,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[13,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[13,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[13,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[13,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[13,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[13,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[13,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/domp.jpg')

# POM Carbon
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Particulate Organic Carbon [mg $C$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[14,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[14,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[14,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[14,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[14,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[14,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[14,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[14,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[14,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[14,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[14,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[14,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[14,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[14,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/pomc.jpg')

# POM Nitrogen
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Particulate Organic Nitrogen [mmol $N$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[15,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[15,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[15,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[15,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[15,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[15,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[15,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[15,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[15,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[15,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[15,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[15,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[15,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[15,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/pomn.jpg')

# POM Phosphorus
fig,ax = plt.subplots(2,4,figsize=[16,10], sharex=True)
handles,labels = ax[0,0].get_legend_handles_labels()
plt.suptitle('Particulate Organic Phosphorus [mmol $P$ $m^{-3}$]')

ax[0,0].set_title('Depth - 0m')
ax[0,0].plot(days, bfm17_daily[16,0,:], '-b', label='GLOBE')
ax[0,0].plot(days, data_fortran[16,0,:], '-.k', label='Fortran')

ax[0,1].set_title('Depth - 25m')
ax[0,1].plot(days, bfm17_daily[16,24,:], '-b', label='GLOBE')
ax[0,1].plot(days, data_fortran[16,24,:], '-.k', label='Fortran')

ax[0,2].set_title('Depth - 50m')
ax[0,2].plot(days, bfm17_daily[16,49,:], '-b', label='GLOBE')
ax[0,2].plot(days, data_fortran[16,49,:], '-.k', label='Fortran')

ax[0,3].set_title('Depth - 75m')
ax[0,3].plot(days, bfm17_daily[16,74,:], '-b', label='GLOBE')
ax[0,3].plot(days, data_fortran[16,74,:], '-.k', label='Fortran')

ax[1,0].set_title('Depth - 100m')
ax[1,0].plot(days, bfm17_daily[16,99,:], '-b', label='GLOBE')
ax[1,0].plot(days, data_fortran[16,99,:], '-.k', label='Fortran')

ax[1,1].set_title('Depth - 125m')
ax[1,1].plot(days, bfm17_daily[16,124,:], '-b', label='GLOBE')
ax[1,1].plot(days, data_fortran[16,124,:], '-.k', label='Fortran')

ax[1,2].set_title('Depth - 150m')
ax[1,2].plot(days, bfm17_daily[16,149,:], '-b', label='GLOBE')
ax[1,2].plot(days, data_fortran[16,149,:], '-.k', label='Fortran')

ax[1,3].remove()
plt.tight_layout()
plt.savefig(os.getcwd() + '/tests/bfm17/figures/pomp.jpg')


# ----------------------------------------------------------------------------------------------------
# 1D Field Plots
# ----------------------------------------------------------------------------------------------------
# rmse_data, nrmse_data = nrmse(avg_data_fortran,bfm17_monthly[:,:,:])    # 2nd year
rmse_data, nrmse_data = nrmse(avg_data_fortran,bfm17_monthly)    # 2nd year
fields = ['o2', 'no3', 'nh4', 'po4', 'phyto_c', 'phyto_n', 'phyto_p', 'phyto_l',
          'zoo_c', 'zoo_n', 'zoo_p', 'dom_c', 'dom_n', 'dom_p', 'pom_c', 'pom_n', 'pom_p']

print('-------------------------------------------------')
print('NRMSE (%)')
print('-------------------------------------------------')
for i in range(0,len(fields)):
    print(fields[i], '--', nrmse_data[i])
print('-------------------------------------------------')
print('RMSE')
print('-------------------------------------------------')
for i in range(0,len(fields)):
    print(fields[i], '--', rmse_data[i])
print('-------------------------------------------------')

# Nutrients
fig,axes = plt.subplots(2,4, figsize=[16,10])
plt.subplot(2,4,1)
plt.imshow(bfm17_monthly[0,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150])
plt.clim(180,230)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('GLOBE - Oxygen')

plt.subplot(2,4,2)
plt.imshow(bfm17_monthly[1,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,2.5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Nitrate')

plt.subplot(2,4,3)
plt.imshow(bfm17_monthly[3,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.075)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Phosphate')

plt.subplot(2,4,4)
plt.imshow(bfm17_monthly[2,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,2.5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Ammonium')

plt.subplot(2,4,5)
plt.imshow(avg_data_fortran[0,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150])
plt.clim(180,230)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('BFM17 - Oxygen')

plt.subplot(2,4,6)
plt.imshow(avg_data_fortran[1,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,2.5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Nitrate')

plt.subplot(2,4,7)
plt.imshow(avg_data_fortran[3,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.075)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Phosphate')

plt.subplot(2,4,8)
plt.imshow(avg_data_fortran[2,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,2.5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Ammonium')

plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/nutrients-1d.jpg')


# Phytoplankton
fig,axes = plt.subplots(2,4, figsize=[16,10])
plt.subplot(2,4,1)
plt.imshow(bfm17_monthly[4,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150])
plt.clim(0,25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('GLOBE - Phyto Carbon')

plt.subplot(2,4,2)
plt.imshow(bfm17_monthly[5,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Phyto Nitrogen')

plt.subplot(2,4,3)
plt.imshow(bfm17_monthly[6,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.02)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Phyto Phosphorus')

plt.subplot(2,4,4)
plt.imshow(bfm17_monthly[7,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Phyto Chlorophyll-a')

plt.subplot(2,4,5)
plt.imshow(avg_data_fortran[4,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150])
plt.clim(0,25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('BFM17 - Phyto Carbon')

plt.subplot(2,4,6)
plt.imshow(avg_data_fortran[5,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Phyto Nitrogen')

plt.subplot(2,4,7)
plt.imshow(avg_data_fortran[6,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.02)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Phyto Phosphorus')

plt.subplot(2,4,8)
plt.imshow(avg_data_fortran[7,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Phyto Chlorophyll-a')

plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phyto-1d.jpg')

# phyto_data = [bfm17_monthly[4,:,:],  bfm17_monthly[5,:,:],  bfm17_monthly[6,:,:],  bfm17_monthly[7,:,:],
#               avg_data_fortran[4, :, :], avg_data_fortran[5, :, :], avg_data_fortran[6, :, :], avg_data_fortran[7, :, :]]

# phyto_titles = ['GLOBE - Phyto Carbon', 'GLOBE - Phyto Nitrogen', 'GLOBE - Phyto Phosphorus', 'GLOBE - Phyto Chlorophyll-a',
#                 'BFM17 - Phyto Carbon', 'BFM17 - Phyto Nitrogen', 'BFM17 - Phyto Phosphorus', 'BFM17 - Phyto Chlorophyll-a']

# limits = [[0,25], [0,0.25], [0,0.02], [0,0.25],
#           [0,25], [0,0.25], [0,0.02], [0,0.25]]

# xticks = [0.5, 2.5, 4.5, 6.5, 8.5, 10.5]

# fig,axes = plt.subplots(2,4, figsize=[16,10])
# for i,ax in enumerate(axes.flat):
#     im = ax.imshow(phyto_data[i], extent=[0,12,150,0], aspect='auto',cmap=cmaps.plasma, vmin=limits[i][0], vmax=limits[i][1])
#     ax.set_xticks(xticks)
#     if i < 4:
#         ax.set_xticklabels(['', '', '', '', '', ''])
#     else:
#         ax.set_xticklabels(['J', 'M', 'M', 'J', 'S', 'N'])
#         ax.set_xlabel('Month')

#     ax.set_yticks([0, 50, 100, 150])

#     if i % 4 != 0:  ax.set_yticklabels(['', '', '', ''])
#     else:   ax.set_ylabel('Depth (m)')

#     ax.set_title(phyto_titles[i])

#     # divider = make_axes_locatable(ax)
#     # cax = divider.append_axes("right", size="5%", pad=0.025)

#     # fig.colorbar(im, cax=cax)

plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/phyto-1d.jpg')
plt.close()

# Zooplankton
fig,axes = plt.subplots(2,3, figsize=[16,10])
plt.subplot(2,3,1)
plt.imshow(bfm17_monthly[8,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150])
plt.clim(0,25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('GLOBE - Zoo Carbon')

plt.subplot(2,3,2)
plt.imshow(bfm17_monthly[9,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Zoo Nitrogen')

plt.subplot(2,3,3)
plt.imshow(bfm17_monthly[10,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.015)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - Zoo Phosphorus')

plt.subplot(2,3,4)
plt.imshow(avg_data_fortran[8,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150])
plt.clim(0,25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('BFM17 - Zoo Carbon')

plt.subplot(2,3,5)
plt.imshow(avg_data_fortran[9,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.25)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Zoo Nitrogen')

plt.subplot(2,3,6)
plt.imshow(avg_data_fortran[10,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.015)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - Zoo Phosphorus')
plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/zoo-1d.jpg')


# DOM
fig,axes = plt.subplots(2,3, figsize=[16,10])
plt.subplot(2,3,1)
plt.imshow(bfm17_monthly[11,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150])
plt.clim(0,200)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('GLOBE - DOM Carbon')

plt.subplot(2,3,2)
plt.imshow(bfm17_monthly[12,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.075)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - DOM Nitrogen')

plt.subplot(2,3,3)
plt.imshow(bfm17_monthly[13,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.0035)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - DOM Phosphorus')

plt.subplot(2,3,4)
plt.imshow(avg_data_fortran[11,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150])
plt.clim(0,200)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('BFM17 - DOM Carbon')

plt.subplot(2,3,5)
plt.imshow(avg_data_fortran[12,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.075)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - DOM Nitrogen')

plt.subplot(2,3,6)
plt.imshow(avg_data_fortran[13,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.0035)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - DOM Phosphorus')
plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/dom-1d.jpg')


# POM
fig,axes = plt.subplots(2,3, figsize=[16,10])
plt.subplot(2,3,1)
plt.imshow(bfm17_monthly[14,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150])
plt.clim(0,1.75)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('GLOBE - POM Carbon')

plt.subplot(2,3,2)
plt.imshow(bfm17_monthly[15,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.025)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - POM Nitrogen')

plt.subplot(2,3,3)
plt.imshow(bfm17_monthly[16,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['','','','','',''])
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.001)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('GLOBE - POM Phosphorus')

plt.subplot(2,3,4)
plt.imshow(avg_data_fortran[14,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,1.75)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_ylabel('Depth (m)')
ax.set_title('BFM17 - POM Carbon')

plt.subplot(2,3,5)
plt.imshow(avg_data_fortran[15,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.025)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - POM Nitrogen')

plt.subplot(2,3,6)
plt.imshow(avg_data_fortran[16,:,:],extent=[0,12,150,0],aspect='auto',cmap=cmaps.plasma)
ax = plt.gca()
plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
plt.xlabel('Month')
plt.yticks([0,50,100,150],['','','',''])
plt.clim(0,0.001)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.025)
plt.colorbar(cax=cax)
ax.set_title('BFM17 - POM Phosphorus')
plt.tight_layout(h_pad=0.75, w_pad=0.75)
plt.savefig(os.getcwd() + '/tests/bfm17/figures/pom-1d.jpg')

