from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os


path = os.getcwd() + '/tests/cobalt/data/ocean.stats.nc'
variables = nc.Dataset(path)
variables = variables.variables
x=1