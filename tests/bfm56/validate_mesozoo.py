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
    std = np.zeros(len(check))
    for i in range(0,len(check)):
        avg[i] = np.abs(np.mean(check[i,:,:]))
        dif[i] = np.max(check[i,:,:]) - np.min(check[i,:,:])
        rms[i] = np.power( np.mean( np.power( check[i,:,:]-comp[i,:,:], 2 ) )   ,0.5)
        max[i] = np.max(check[i,:,:])
        std[i] = np.std(check[i,:,:])
        # nrmse = 100*rms/avg    
        # nrmse = 100*rms/max
    nrmse = 100*rms/(std + 1.E-20)

    return rms, nrmse

# ----------------------------------------------------------------------------------------------------
# Extract model results
# ----------------------------------------------------------------------------------------------------
# GLOBE
path = os.getcwd() + '/tests/bfm56/data/concentration_bfm56_arbitrary_mesozoo2.npz'
globe_bfm = np.load(path, allow_pickle=True)

globe_daily = globe_bfm["daily"]
globe_monthly = globe_bfm["monthly"]

# Tracer Indices
path = os.getcwd() + '/tests/bfm56/data/tracer_indices_bfm56_arbitrary_mesozoo2.npz'
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
path = os.getcwd() + '/tests/bfm56/data/bfm56_pom1d_arbitrary_mesozoo2.nc'
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

# bfm56_monthly = np.zeros((50,150,120))
# for spec in range(0,50):
#     for year in range(0,10):
#         for month in range(0,12):
#             for day in range(0,30):
#                 bfm56_monthly[spec,:,month + (year*12)] = bfm56_monthly[spec,:,month + (year*12)] + bfm56_daily[spec,:,(day + (month*30) + (year*360))]
# bfm56_monthly = bfm56_monthly/30

bfm56_daily = bfm56_daily[:,:,:1800]
# bfm56_monthly = bfm56_monthly[:,:,:60]
days = np.linspace(0,29,30)
xticks = [0,10,20,30]

fig,ax = plt.subplots(4,3,figsize=[15,10], sharex=True)

ax[0,0].set_title('(a)')
ax[0,0].plot(days, globe_daily[mesoz2c,0,:30], '-b', label='GLOBE')
ax[0,0].plot(days, z4c[0,:30], '-.k', label='BFM56')
ax[0,0].set_xlim([0,30])
ax[0,0].set_ylabel("mg C ${m^{-3}}$")

ax[1,0].set_title('(d)')
ax[1,0].plot(days, globe_daily[mesoz2c,49,:30], '-b', label='GLOBE')
ax[1,0].plot(days, z4c[49,:30], '-.k', label='BFM56')
ax[1,0].set_xlim([0,30])
ax[1,0].set_ylabel("mg C ${m^{-3}}$")

ax[2,0].set_title('(g)')
ax[2,0].plot(days, globe_daily[mesoz2c,99,:30], '-b', label='GLOBE')
ax[2,0].plot(days, z4c[99,:30], '-.k', label='BFM56')
ax[2,0].set_xlim([0,30])
ax[2,0].set_ylabel("mg C ${m^{-3}}$")

ax[3,0].set_title('(j)')
ax[3,0].plot(days, globe_daily[mesoz2c,149,:30], '-b', label='GLOBE')
ax[3,0].plot(days, z4c[149,:30], '-.k', label='BFM56')
ax[3,0].set_xticks(xticks)
ax[3,0].set_xlim([0,30])
ax[3,0].set_xlabel('Time [day]')
ax[3,0].set_ylabel("mg C ${m^{-3}}$")

ax[0,1].set_title('(b)')
ax[0,1].plot(days, globe_daily[mesoz2n,0,:30], '-b', label='GLOBE')
ax[0,1].plot(days, z4n[0,:30], '-.k', label='BFM56')
ax[0,1].set_xlim([0,30])
ax[0,1].set_ylabel("mmol N ${m^{-3}}$")

ax[1,1].set_title('(e)')
ax[1,1].plot(days, globe_daily[mesoz2n,49,:30], '-b', label='GLOBE')
ax[1,1].plot(days, z4n[49,:30], '-.k', label='BFM56')
ax[1,1].set_xlim([0,30])
ax[1,1].set_ylabel("mmol N ${m^{-3}}$")

ax[2,1].set_title('(h)')
ax[2,1].plot(days, globe_daily[mesoz2n,99,:30], '-b', label='GLOBE')
ax[2,1].plot(days, z4n[99,:30], '-.k', label='BFM56')
ax[2,1].set_xlim([0,30])
ax[2,1].set_ylabel("mmol N ${m^{-3}}$")

ax[3,1].set_title('(k)')
ax[3,1].plot(days, globe_daily[mesoz2n,149,:30], '-b', label='GLOBE')
ax[3,1].plot(days, z4n[149,:30], '-.k', label='BFM56')
ax[3,1].set_xticks(xticks)
ax[3,1].set_xlim([0,30])
ax[3,1].set_xlabel('Time [day]')
ax[3,1].set_ylabel("mmol N ${m^{-3}}$")

ax[0,2].set_title('(c)')
ax[0,2].plot(days, globe_daily[mesoz2p,0,:30], '-b', label='GLOBE')
ax[0,2].plot(days, z4p[0,:30], '-.k', label='BFM56')
ax[0,2].set_xlim([0,30])
ax[0,2].set_ylabel("mmol P ${m^{-3}}$")

ax[1,2].set_title('(f)')
ax[1,2].plot(days, globe_daily[mesoz2p,49,:30], '-b', label='GLOBE')
ax[1,2].plot(days, z4p[49,:30], '-.k', label='BFM56')
ax[1,2].set_xlim([0,30])
ax[1,2].set_ylabel("mmol P ${m^{-3}}$")

ax[2,2].set_title('(i)')
ax[2,2].plot(days, globe_daily[mesoz2p,99,:30], '-b', label='GLOBE')
ax[2,2].plot(days, z4p[99,:30], '-.k', label='BFM56')
ax[2,2].set_xlim([0,30])
ax[2,2].set_ylabel("mmol P ${m^{-3}}$")

ax[3,2].set_title('(l)')
ax[3,2].plot(days, globe_daily[mesoz2p,149,:30], '-b', label='GLOBE')
ax[3,2].plot(days, z4p[149,:30], '-.k', label='BFM56')
ax[3,2].set_xticks(xticks)
ax[3,2].set_xlim([0,30])
ax[3,2].set_xlabel('Time [day]')
ax[3,2].set_ylabel("mmol P ${m^{-3}}$")

handles,labels = ax[0,0].get_legend_handles_labels()
# fig.legend(handles,labels, loc='lower center', ncol=2)

fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
plt.savefig(os.getcwd() + '/tests/bfm56/figures/validate_mesozoo.jpg')
plt.close()

rmse_data, nrmse_data = nrmse(bfm56_daily[:,:30],globe_daily[:,:30])    # 2nd year
print('-------------------------------------------------')
print('NRMSE (%)')
print('-------------------------------------------------')
for i in range(0,len(globe_daily)):
    print(nrmse_data[i])
print('-------------------------------------------------')
print('RMSE')
print('-------------------------------------------------')
for i in range(0,len(globe_daily)):
    print(rmse_data[i])
print('-------------------------------------------------')