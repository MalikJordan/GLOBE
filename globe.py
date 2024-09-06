import os
import sys
import numpy as np
from setup.initialize import import_model
from functions.other_functions import *
from functions import rates
from functions.bgc_rate_eqns import bgc_rate_eqns
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------
# Import and initialize model
# ----------------------------------------------------------------------------------------------------
# check_file = False
# first_check = True
# while not check_file:
#     if first_check:
#         file = input("\nEnter the name of the yaml input file. Press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     else:
#         file = input("\nInput file not found. Re-enter the name of the yaml input file or press 'Q' to quit. \nEx: model_description.yaml \n\n")
#     if file == 'Q' or file == 'q':
#         sys.exit()
#     check_file = os.path.exists(file)
#     first_check = False
file = 'npzd.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

# ----------------------------------------------------------------------------------------------------
# Begin simulation
# ----------------------------------------------------------------------------------------------------
end_time = 86400 * parameters["simulation"]["num_days"]
time = np.arange(0, end_time + parameters["simulation"]["timestep"], parameters["simulation"]["timestep"])

for i in range(1,len(time)):
    bgc_rate_eqns(time[i-1], file_path, tracers)

#     conc += bgc_rates * parameters["simulation"]["timestep"]
#     solution[:,i] = conc

