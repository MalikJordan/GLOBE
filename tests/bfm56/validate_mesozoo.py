from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os
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

mesoz2c = globe_daily[tracer_indices["mesozoo2"][0]]
mesoz2n = globe_daily[tracer_indices["mesozoo2"][1]]
mesoz2p = globe_daily[tracer_indices["mesozoo2"][2]]

# pyPOM1D-BFM
path = os.getcwd() + '/tests/bfm56/data/bfm56_pom1d_arbitrary_mesozoo2.nc'
variables = nc.Dataset(path)
variables = variables.variables

z4c = np.asarray(variables['Z4c'][:]).transpose()   # omnivorous mesozooplankton
z4n = np.asarray(variables['Z4n'][:]).transpose()
z4p = np.asarray(variables['Z4p'][:]).transpose()

days = np.linspace(0,29,30)
xticks = [0,10,20,30]

fig,ax = plt.subplots(4,3,figsize=[15,10], sharex=True)

ax[0,0].set_title('(a)')
ax[0,0].plot(days, mesoz2c[0,:30], '-b', label='GLOBE')
ax[0,0].plot(days, z4c[0,:30], '-.k', label='BFM56')
ax[0,0].set_ylabel("mg C ${m^{-3}}$")

ax[1,0].set_title('(d)')
ax[1,0].plot(days, mesoz2c[49,:30], '-b', label='GLOBE')
ax[1,0].plot(days, z4c[49,:30], '-.k', label='BFM56')
ax[1,0].set_ylabel("mg C ${m^{-3}}$")

ax[2,0].set_title('(g)')
ax[2,0].plot(days, mesoz2c[99,:30], '-b', label='GLOBE')
ax[2,0].plot(days, z4c[99,:30], '-.k', label='BFM56')
ax[2,0].set_ylabel("mg C ${m^{-3}}$")

ax[3,0].set_title('(j)')
ax[3,0].plot(days, mesoz2c[149,:30], '-b', label='GLOBE')
ax[3,0].plot(days, z4c[149,:30], '-.k', label='BFM56')
ax[3,0].set_xticks(xticks)
ax[3,0].set_ylabel("mg C ${m^{-3}}$")

ax[0,1].set_title('(b)')
ax[0,1].plot(days, mesoz2n[0,:30], '-b', label='GLOBE')
ax[0,1].plot(days, z4n[0,:30], '-.k', label='BFM56')
ax[0,1].set_ylabel("mmol N ${m^{-3}}$")

ax[1,1].set_title('(e)')
ax[1,1].plot(days, mesoz2n[49,:30], '-b', label='GLOBE')
ax[1,1].plot(days, z4n[49,:30], '-.k', label='BFM56')
ax[1,1].set_ylabel("mmol N ${m^{-3}}$")

ax[2,1].set_title('(h)')
ax[2,1].plot(days, mesoz2n[99,:30], '-b', label='GLOBE')
ax[2,1].plot(days, z4n[99,:30], '-.k', label='BFM56')
ax[2,1].set_ylabel("mmol N ${m^{-3}}$")

ax[3,1].set_title('(k)')
ax[3,1].plot(days, mesoz2n[149,:30], '-b', label='GLOBE')
ax[3,1].plot(days, z4n[149,:30], '-.k', label='BFM56')
ax[3,1].set_xticks(xticks)
ax[3,1].set_ylabel("mmol N ${m^{-3}}$")

ax[0,2].set_title('(c)')
ax[0,2].plot(days, mesoz2p[0,:30], '-b', label='GLOBE')
ax[0,2].plot(days, z4p[0,:30], '-.k', label='BFM56')
ax[0,2].set_ylabel("mmol P ${m^{-3}}$")

ax[1,2].set_title('(f)')
ax[1,2].plot(days, mesoz2p[49,:30], '-b', label='GLOBE')
ax[1,2].plot(days, z4p[49,:30], '-.k', label='BFM56')
ax[3,2].set_ylabel("mmol P ${m^{-3}}$")

ax[2,2].set_title('(i)')
ax[2,2].plot(days, mesoz2p[99,:30], '-b', label='GLOBE')
ax[2,2].plot(days, z4p[99,:30], '-.k', label='BFM56')
ax[3,2].set_ylabel("mmol P ${m^{-3}}$")

ax[3,2].set_title('(l)')
ax[3,2].plot(days, mesoz2p[149,:30], '-b', label='GLOBE')
ax[3,2].plot(days, z4p[149,:30], '-.k', label='BFM56')
ax[3,2].set_xticks(xticks)
ax[3,2].set_ylabel("mmol P ${m^{-3}}$")

handles,labels = ax[0,0].get_legend_handles_labels()
fig.legend(handles,labels, loc='lower center', ncol=2)

fig.tight_layout(h_pad=2.5,w_pad=2.5,rect=[0,0.025,1,1])
plt.savefig(os.getcwd() + '/tests/bfm56/figures/validate_mesozoo.jpg')
plt.close()