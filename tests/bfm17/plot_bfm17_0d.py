import netCDF4 as nc
import os
import numpy as np
from matplotlib import pyplot as plt

def nrmse(check,comp):
    
    avg = np.zeros(len(check))
    dif = np.zeros(len(check))
    rms = np.zeros(len(check))
    max = np.zeros(len(check))
    std = np.zeros(len(check))
    for i in range(0,len(check)):
        avg[i] = np.abs(np.mean(check[i,:]))
        dif[i] = np.max(check[i,:]) - np.min(check[i,:])
        rms[i] = np.power( np.mean( np.power( check[i,:]-comp[i,:], 2 ) )   ,0.5)
        max[i] = np.max(check[i,:])
        std[i] = np.std(check[i,:])
    # nrmse = 100*rms/avg    
    # nrmse = 100*rms/max
    nrmse = 100*rms/(std + 1.E-20)

    return rms, nrmse

folder = os.getcwd() + '/tests/bfm17'

# Get Fortran Data --------------------------------------------------------------------
path = os.getcwd() + '/tests/bfm17/data/BFM_standalone_pelagic.nc'
variables = nc.Dataset(path).variables

time = np.asarray(variables['time'])
o2_bfm17 = np.asarray(variables['O2o'][:360])
no3_bfm17 = np.asarray(variables['N3n'][:360])
nh4_bfm17 = np.asarray(variables['N4n'][:360])
po4_bfm17 = np.asarray(variables['N1p'][:360])
hs_bfm17 = np.asarray(variables['N6r'][:360])
pc_bfm17 = np.asarray(variables['P2c'][:360])
pn_bfm17 = np.asarray(variables['P2n'][:360])
pp_bfm17 = np.asarray(variables['P2p'][:360])
pl_bfm17 = np.asarray(variables['P2l'][:360])
zc_bfm17 = np.asarray(variables['Z5c'][:360])
zn_bfm17 = np.asarray(variables['Z5n'][:360])
zp_bfm17 = np.asarray(variables['Z5p'][:360])
domc_bfm17 = np.asarray(variables['R1c'][:360])
domn_bfm17 = np.asarray(variables['R1n'][:360])
domp_bfm17 = np.asarray(variables['R1p'][:360])
pomc_bfm17 = np.asarray(variables['R6c'][:360])
pomn_bfm17 = np.asarray(variables['R6n'][:360])
pomp_bfm17 = np.asarray(variables['R6p'][:360])

bfm17 = np.zeros((18,360))
bfm17[0,:] = o2_bfm17[:,0]
bfm17[1,:] = no3_bfm17[:,0]
bfm17[2,:] = nh4_bfm17[:,0]
bfm17[3,:] = po4_bfm17[:,0]
bfm17[4,:] = hs_bfm17[:,0]
bfm17[5,:] = pc_bfm17[:,0]
bfm17[6,:] = pn_bfm17[:,0]
bfm17[7,:] = pp_bfm17[:,0]
bfm17[8,:] = pl_bfm17[:,0]
bfm17[9,:] = zc_bfm17[:,0]
bfm17[10,:] = zn_bfm17[:,0]
bfm17[11,:] = zp_bfm17[:,0]
bfm17[12,:] = domc_bfm17[:,0]
bfm17[13,:] = domn_bfm17[:,0]
bfm17[14,:] = domp_bfm17[:,0]
bfm17[15,:] = pomc_bfm17[:,0]
bfm17[16,:] = pomn_bfm17[:,0]
bfm17[17,:] = pomp_bfm17[:,0]


# Get GLOBE Data --------------------------------------------------------------------
# Load solution
# path = os.getcwd() + "/tests/bfm17/data/concentration_bfm17_0d.npz"
path = os.getcwd() + "/concentration_bfm17_0d.npz"
model = np.load(path, allow_pickle=True)
globe_daily = model["daily"]
globe_monthly = model["monthly"]

globe_daily = globe_daily[...,:360]
globe_monthly = globe_monthly[...,:12]

# Load tracer indices
# path = os.getcwd() + "/tests/bfm17/data/tracer_indices_bfm17_0d.npz"
path = os.getcwd() + "/tracer_indices_bfm17_0d.npz"
indices = np.load(path)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

o2 = tracer_indices["o2"][0]
no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
po4 = tracer_indices["po4"][0]
hs = tracer_indices["hs"][0]
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
xlabel = ['J','M','M','J','S','N']
xticks = [15,75,135,195,255,315]
days = np.linspace(0,359,360)
# Create Plots --------------------------------------------------------------------
fig,axs = plt.subplots(2,4,figsize=(15,8), sharex=True)
axs[0,0].plot(days,globe_daily[o2,0],'-b',label='GLOBE')
axs[0,0].plot(days,o2_bfm17,'-.k',label='BFM17')
# axs[0,0].set_title("(a) Oxygen")
axs[0,0].set_title("(a)")
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mmol O ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[no3,0],'-b',label='GLOBE')
axs[0,1].plot(days,no3_bfm17,'-.k',label='BFM17')
# axs[0,1].set_title("(b) Nitrate")
axs[0,1].set_title("(b)")
axs[0,1].set_xlim([0,360])
axs[0,1].set_ylabel("mmol N ${m^{-3}}$")

axs[0,2].plot(days,globe_daily[nh4,0],'-b',label='GLOBE')
axs[0,2].plot(days,nh4_bfm17,'-.k',label='BFM17')
# axs[0,2].set_title("(c) Ammonium")
axs[0,2].set_title("(c)")
axs[0,2].set_xlim([0,360])
axs[0,2].set_ylabel("mmol N ${m^{-3}}$")

axs[0,3].plot(days,globe_daily[po4,0],'-b',label='GLOBE')
axs[0,3].plot(days,po4_bfm17,'-.k',label='BFM17')
# axs[0,3].set_title("(d) Phosphate")
axs[0,3].set_title("(d)")
axs[0,3].set_xlim([0,360])
axs[0,3].set_ylabel("mmol P ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[pc,0],'-b',label='GLOBE')
axs[1,0].plot(days,pc_bfm17,'-.k',label='BFM17')
# axs[1,0].set_title("(e) Phytoplankton")
axs[1,0].set_title("(e)")
axs[1,0].set_xlabel("Time [month]")
axs[1,0].set_xticks(xticks,xlabel)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mg C ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[zc,0],'-b',label='GLOBE')
axs[1,1].plot(days,zc_bfm17,'-.k',label='BFM17')
# axs[1,1].set_title("(f) Zooplankton")
axs[1,1].set_title("(f)")
axs[1,1].set_xlabel("Time [month]")
axs[1,1].set_xticks(xticks,xlabel)
axs[1,1].set_xlim([0,360])
axs[1,1].set_ylabel("mg C ${m^{-3}}$")

axs[1,2].plot(days,globe_daily[domc,0],'-b',label='GLOBE')
axs[1,2].plot(days,domc_bfm17,'-.k',label='BFM17')
# axs[1,2].set_title("(g) Dissolved Organic Carbon")
axs[1,2].set_title("(g)")
axs[1,2].set_xlabel("Time [month]")
axs[1,2].set_xticks(xticks,xlabel)
axs[1,2].set_xlim([0,360])
axs[1,2].set_ylabel("mg C ${m^{-3}}$")

axs[1,3].plot(days,globe_daily[pomc,0],'-b',label='GLOBE')
axs[1,3].plot(days,pomc_bfm17,'-.k',label='BFM17')
# axs[1,3].set_title("(h) Particulate Organic Carbon")
axs[1,3].set_title("(h)")
axs[1,3].set_xlabel("Time [month]")
axs[1,3].set_xticks(xticks,xlabel)
axs[1,3].set_xlim([0,360])
axs[1,3].set_ylabel("mg C ${m^{-3}}$")

handles, labels = axs[0,0].get_legend_handles_labels()
fig.legend(handles,labels, loc='lower center', ncol=2)

# fig.suptitle("Nutrients")
fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
loc = os.path.join(folder + "/figures","bfm17.jpg")
plt.savefig(loc)

# Nutrients
fig, axs = plt.subplots(2,2,figsize=(10,10),sharex=True)

axs[0,0].plot(days,globe_daily[o2,0],'-b',label='GLOBE')
axs[0,0].plot(days,o2_bfm17,'-.k',label='BFM17')
axs[0,0].set_title("(a) Oxygen")
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mmol O ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[no3,0],'-b',label='GLOBE')
axs[0,1].plot(days,no3_bfm17,'-.k',label='BFM17')
axs[0,1].set_title("(b) Nitrate")
axs[0,1].set_xlim([0,360])
axs[0,1].set_ylabel("mmol N ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[nh4,0],'-b',label='GLOBE')
axs[1,0].plot(days,nh4_bfm17,'-.k',label='BFM17')
axs[1,0].set_title("(c) Ammonium")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mmol N ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[po4,0],'-b',label='GLOBE')
axs[1,1].plot(days,po4_bfm17,'-.k',label='BFM17')
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

axs[0,0].plot(days,globe_daily[pc,0],'-b',label='GLOBE')
axs[0,0].plot(days,pc_bfm17,'-.k',label='BFM17')
axs[0,0].set_title("(a) Phytoplankton -- Carbon")
axs[0,0].set_xlim([0,360])
axs[0,0].set_ylabel("mg C ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[zc,0],'-b',label='GLOBE')
axs[0,1].plot(days,zc_bfm17,'-.k',label='BFM17')
axs[0,1].set_title("(b) Zooplankton -- Carbon")
axs[0,1].set_xlim([0,360])
axs[0,1].set_ylabel("mg C ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[domc,0],'-b',label='GLOBE')
axs[1,0].plot(days,domc_bfm17,'-.k',label='BFM17')
axs[1,0].set_title("(c) Dissolved Organic Carbon")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
axs[1,0].set_xlim([0,360])
axs[1,0].set_ylabel("mg C ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[pomc,0],'-b',label='GLOBE')
axs[1,1].plot(days,pomc_bfm17,'-.k',label='BFM17')
axs[1,1].set_title("(d) Particulate Organic Carbon")
axs[1,1].set_xlabel("Time [months]")
axs[1,1].set_xticks(xticks,xlabel)
axs[1,1].set_xlim([0,360])
axs[1,1].set_ylabel("mg C ${m^{-3}}$")

handles, labels = axs[0,0].get_legend_handles_labels()
# fig.legend(handles,labels, loc='lower center', ncol=2)

# fig.suptitle("Organic Carbon")
fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
org = os.path.join(folder + "/figures","organic.jpg")
plt.savefig(org)

error = np.zeros((8,360))
# error[0] = 100 * (o2_bfm17[:360,0] - globe_daily[o2,0,:360]) / (o2_bfm17[:360,0] + 1.E-20)
# error[1] = 100 * (po4_bfm17[:360,0] - globe_daily[po4,0,:360]) / (po4_bfm17[:360,0] + 1.E-20)
# error[2] = 100 * (no3_bfm17[:360,0] - globe_daily[no3,0,:360]) / (no3_bfm17[:360,0] + 1.E-20)
# error[3] = 100 * (nh4_bfm17[:360,0] - globe_daily[nh4,0,:360]) / (nh4_bfm17[:360,0] + 1.E-20)
# error[4] = 100 * (pc_bfm17[:360,0] - globe_daily[pc,0,:360]) / (pc_bfm17[:360,0] + 1.E-20)
# error[5] = 100 * (zc_bfm17[:360,0] - globe_daily[zc,0,:360]) / (zc_bfm17[:360,0] + 1.E-20)
# error[6] = 100 * (domc_bfm17[:360,0] - globe_daily[domc,0,:360]) / (domc_bfm17[:360,0] + 1.E-20)
# error[7] = 100 * (pomc_bfm17[:360,0] - globe_daily[pomc,0,:360]) / (pomc_bfm17[:360,0] + 1.E-20)

for i in range(360):
    error[0,i] = 100 * (o2_bfm17[i,0] - globe_daily[o2,0,i]) / (o2_bfm17[i,0] + 1.E-20)
    error[1,i] = 100 * (po4_bfm17[i,0] - globe_daily[po4,0,i]) / (po4_bfm17[i,0] + 1.E-20)
    error[2,i] = 100 * (no3_bfm17[i,0] - globe_daily[no3,0,i]) / (no3_bfm17[i,0] + 1.E-20)
    error[3,i] = 100 * (nh4_bfm17[i,0] - globe_daily[nh4,0,i]) / (nh4_bfm17[i,0] + 1.E-20)
    error[4,i] = 100 * (pc_bfm17[i,0] - globe_daily[pc,0,i]) / (pc_bfm17[i,0] + 1.E-20)
    error[5,i] = 100 * (zc_bfm17[i,0] - globe_daily[zc,0,i]) / (zc_bfm17[i,0] + 1.E-20)
    error[6,i] = 100 * (domc_bfm17[i,0] - globe_daily[domc,0,i]) / (domc_bfm17[i,0] + 1.E-20)
    error[7,i] = 100 * (pomc_bfm17[i,0] - globe_daily[pomc,0,i]) / (pomc_bfm17[i,0] + 1.E-20)
print(np.max(error))

check = np.zeros((8,360))
comp = np.zeros((8,360))

check[0,:] = o2_bfm17[:,0]
check[1,:] = po4_bfm17[:,0]
check[2,:] = no3_bfm17[:,0]
check[3,:] = nh4_bfm17[:,0]
check[4,:] = pc_bfm17[:,0]
check[5,:] = zc_bfm17[:,0]
check[6,:] = domc_bfm17[:,0]
check[7,:] = pomc_bfm17[:,0]

comp[0,:] = globe_daily[o2,0,:]
comp[1,:] = globe_daily[po4,0,:]
comp[2,:] = globe_daily[no3,0,:]
comp[3,:] = globe_daily[nh4,0,:]
comp[4,:] = globe_daily[pc,0,:]
comp[5,:] = globe_daily[zc,0,:]
comp[6,:] = globe_daily[domc,0,:]
comp[7,:] = globe_daily[pomc,0,:]

check = np.zeros((18,360))

# rmse_data, nrmse_data = nrmse(check,comp)
# fields = ['o2', 'po4', 'no3', 'nh4', 'phyto_c', 'zoo_c', 'dom_c', 'pom_c']
rmse_data, nrmse_data = nrmse(bfm17,globe_daily[:,0])
fields = ['o2', 'no3', 'nh4', 'po4', 'hs', 'phyto_c', 'phyto_n', 'phyto_p', 'phyto_l', 'zoo_c', 'zoo_n', 'zoo_p', 'dom_c', 'dom_n', 'dom_p', 'pom_c', 'pom_n', 'pom_p']
print(np.max(globe_daily[pc,0,:]))
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
