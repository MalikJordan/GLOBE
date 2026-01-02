import netCDF4 as nc
import os
import numpy as np
from matplotlib import pyplot as plt

path = os.getcwd() + '/tests/bfm17/data/BFM_standalone_pelagic.nc'
variables = nc.Dataset(path).variables

o2 = np.asarray(variables['O2o'][:])
po4 = np.asarray(variables['N1p'][:])
no3 = np.asarray(variables['N3n'][:])
nh4 = np.asarray(variables['N4n'][:])

pc = np.asarray(variables['P2c'][:])
pn = np.asarray(variables['P2n'][:])
pp = np.asarray(variables['P2p'][:])
pl = np.asarray(variables['P2l'][:])
exu = np.asarray(variables['exPPYc'][:])
gpp = np.asarray(variables['ruPPYc'][:])
uptn = np.asarray(variables['ruPPYn'][:])
uptp = np.asarray(variables['ruPPYp'][:])

domc = np.asarray(variables['R1c'][:])
domn = np.asarray(variables['R1n'][:])
domp = np.asarray(variables['R1p'][:])

pomc = np.asarray(variables['R6c'][:])
pomn = np.asarray(variables['R6n'][:])
pomp = np.asarray(variables['R6p'][:])

resz = np.asarray(variables['resZOOc'][:])
nspz = np.asarray(variables['ruZOOc'][:])
zc = np.asarray(variables['Z5c'][:])

remzn = np.asarray(variables['remZOOn'][:])
remzp = np.asarray(variables['remZOOp'])

fI = np.asarray(variables['eiPPY_iiP2_'][:])
ETW = np.asarray(variables['ETW'][:])
EIR = np.asarray(variables['EIR'][:])

x = np.linspace(0,729,730)