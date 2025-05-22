import numpy as np
from matplotlib import pyplot as plt
import os

# Folder path
folder = os.getcwd() + '/tests/2npzd-2'

# Load solution
solution = np.load("npzd-0404.npz", allow_pickle=True)
conc = solution["conc"]     # concentratrion matrix
time = solution["time"]     # time array

# Load tracer indices
indices = np.load("tracer_indices-0220.npz")
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

# Convert time array from seconds to months
sec_day = 60 * 60 * 24
time = time/sec_day