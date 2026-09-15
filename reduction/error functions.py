import numpy as np
import copy
from scipy.integrate import solve_ivp
from functions.bgc_rate_eqns import bgc_rate_eqns

def peak(t_min, t_max, solution_full, solution_reduced, tracer_map_full, tracer_map_reduced):

    # Find index associated with t_min and t_max for slicing the list
    for index,time in enumerate(solution_full.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full.t):
        if time > t_max:
            index_t_max_full = index
            break

    for index,time in enumerate(solution_reduced.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced.t):
        if time > t_max:
            index_t_max_reduced = index
            break

    