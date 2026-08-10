import numpy as np
np.set_printoptions(precision=20)

def average(concentration, physical, type):

    timesteps_per_day = 86400/physical["simulation"]["dt"]
    day = 0
    month = 0

    if type == 'concentration':
        daily = np.zeros((concentration.shape[0],concentration.shape[1],physical["simulation"]["days"]))
        monthly = np.zeros((concentration.shape[0],concentration.shape[1],physical["simulation"]["months"]))

        for i in range(0,physical["simulation"]["iters"]-1):
            daily[:,:,day] += concentration[:,:,i] # add to day tally

            if (i != 0) & ((i+1) % timesteps_per_day == 0): # take average at the end of day
                daily[:,:,day] = daily[:,:,day]/timesteps_per_day
                monthly[:,:,month] += daily[:,:,day] # add to month tally

                if (day != 0) & ((day+1) % 30 == 0):
                    monthly[:,:,month] = monthly[:,:,month]/30 # take average at the end of month
                    month += 1 # move to next month

                day += 1 # move to next day

    else:
        daily = np.zeros((concentration.shape[0],physical["simulation"]["days"]))
        monthly = np.zeros((concentration.shape[0],physical["simulation"]["months"]))

        for i in range(0,physical["simulation"]["iters"]-1):
            daily[:,day] += concentration[:,i] # add to day tally

            if (i != 0) & ((i+1) % timesteps_per_day == 0): # take average at the end of day
                daily[:,day] = daily[:,day]/timesteps_per_day
                monthly[:,month] += daily[:,day] # add to month tally

                if (day != 0) & ((day+1) % 30 == 0):
                    monthly[:,month] = monthly[:,month]/30 # take average at the end of month
                    month += 1 # move to next month

                day += 1 # move to next day

    return daily, monthly