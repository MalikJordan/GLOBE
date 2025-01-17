import os
import sys
import numpy as np
from setup.initialize import import_model
from functions.other_functions import *
from functions import rates
from functions.bgc_rate_eqns import bgc_rate_eqns_solveivp
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt

# Initialize Model
file = 'bfm17.yaml'
file_path = os.getcwd() + '/' + file
base_element, parameters, reactions, tracers = import_model(file_path)

# Write tracers as matrix
c0 = []
for t in tracers:
    for i in range(len(tracers[t].composition)):
        c0.append(tracers[t].conc[i,...,0])
c0=np.array(c0,dtype=float)

# Scipy Integration
t_span = [0,86400*365*2]
solution = solve_ivp(bgc_rate_eqns_solveivp, t_span, base_element, parameters, tracers, method='RK23')

x=0