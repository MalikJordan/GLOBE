import netCDF4    
import os
import numpy as np

file = os.getcwd() + '/SimData_yr10.nc'
out = netCDF4.Dataset(file)  # output file
sol_phy = out.variables['Phy'][:]
sol_zoo = out.variables['Zoo'][:]
sol_nut = out.variables['Nut'][:]

phyto = np.asarray(sol_phy[:,0,0,0,0])
zoo = np.asarray(sol_zoo[:,0,0,0,0])
nut = np.asarray(sol_nut[:,0,0,0])

print()