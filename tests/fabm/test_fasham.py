import xarray as xr
import numpy as np
import pandas as pd
import os
import sys
from matplotlib import pyplot as plt
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
from functions.seasonal_cycling import get_temperature, get_salinity, get_sunlight

# # Create envirnmental data file for fasham with constant temperature, salinity, and shortwave radiation
# # Hourly time series for one year
# time = pd.date_range(
#     "1998-01-01 00:00:00",
#     "1999-01-01 00:00:00",
#     freq="1h",
#     inclusive="left"
# )

# df = pd.DataFrame({
#     "date": time.strftime("%Y-%m-%d"),
#     "time": time.strftime("%H:%M:%S"),
#     "shortwave": 100.0,
#     "temperature": 10.0,
#     "salinity": 35.0
# })

# df.to_csv(
#     "env_constant.dat",
#     sep=" ",
#     header=False,
#     index=False
# )

# # Check constant environmental forcings
# env_data_file = os.getcwd() + '/env_constant.dat'
# data = pd.read_csv(
#     env_data_file,
#     sep=r"\s+",
#     header=None,
#     names=["date", "time", "shortwave", "temperature", "salinity"]
# )


# date = data["date"].to_numpy()
# time = data["time"].to_numpy()
# shortwave = data["shortwave"].to_numpy()
# temperature = data["temperature"].to_numpy()
# salinity = data["salinity"].to_numpy()

# x=1

# # Create environmental data file for fasham with seasonal cycling for temperature, salinity, and shortwave radiation
# # Time series for one year with data every 360s
# dt = 360
# num_timesteps = int(365 * 24 * 60 * 60 / 360) + 2
# time_array = np.arange(num_timesteps) * dt

# time = (
#     pd.Timestamp("1998-01-01 00:00:00") + pd.to_timedelta(time_array, unit="s")
# )

# # time = pd.date_range(
# #     "1998-01-01 00:00:00",
# #     "1999-01-01 00:00:00",
# #     freq="360s",
# #     inclusive="left"
# # )

# temperature = np.zeros(len(time))
# salinity = np.zeros(len(time))
# shortwave = np.zeros(len(time))


# for i in range(len(time)):
#     temperature[i] = get_temperature(time_array[i], 7., 15., 0.)
#     salinity[i] = get_salinity(time_array[i], 33.5, 31.)
#     shortwave[i] = get_sunlight(time_array[i], 20., 240., 54.)

# df = pd.DataFrame({
#     "date": time.strftime("%Y-%m-%d"),
#     "time": time.strftime("%H:%M:%S"),
#     "shortwave": shortwave,
#     "temperature": temperature,
#     "salinity": salinity
# })

# df.to_csv(
#     "env_seasonal.dat",
#     sep=" ",
#     header=False,
#     index=False
# )

# x=1








# output_file = os.getcwd() + '/tests/fabm/fasham_data/output_constant.nc'
output_file = os.getcwd() + '/tests/fabm/fasham_data/output_seasonal.nc'
ds = xr.open_dataset(output_file)

# env_data_file = os.getcwd() + '/tests/fabm/fasham_data/env_nns_annual.dat'
# env_data_file = os.getcwd() + '/tests/fabm/fasham_data/env_constant.dat'
env_data_file = os.getcwd() + '/tests/fabm/fasham_data/env_seasonal.dat'
data = pd.read_csv(
    env_data_file,
    sep=r"\s+",
    header=None,
    names=["date", "time", "shortwave", "temperature", "salinity"]
)

date = data["date"].to_numpy()
time = data["time"].to_numpy()
shortwave = data["shortwave"].to_numpy()
temperature = data["temperature"].to_numpy()
salinity = data["salinity"].to_numpy()


fasham = np.zeros((7,87601))
fasham[0] = ds["fasham_nit"].data[:,0,0]
fasham[1] = ds["fasham_amm"].data[:,0,0]
fasham[2] = ds["fasham_bac"].data[:,0,0]
fasham[3] = ds["fasham_phy"].data[:,0,0]
fasham[4] = ds["fasham_zoo"].data[:,0,0]
fasham[5] = ds["fasham_ldn"].data[:,0,0]
fasham[6] = ds["fasham_det"].data[:,0,0]

iters_per_day = 86400/360
fasham_daily = np.zeros((7,366))
day = 0
for i in range(0,87600):
    fasham_daily[:,day] += fasham[:,i]  # add to day tally
    if (i != 0) and ((i+1)%iters_per_day == 0): 
        fasham_daily[:,day] = fasham_daily[:,day]/iters_per_day     # take average at the end of day
        day += 1    # move to next day

fasham_daily = fasham_daily[...,:360]
x = 1

# conc_path = os.getcwd() + '/concentration_fasham_0d_constant.npz'
conc_path = os.getcwd() + '/concentration_fasham_0d_seasonal.npz'
tracer_path = os.getcwd() + '/tracer_indices_fasham_0d.npz'

model = np.load(conc_path, allow_pickle=True)
globe_daily = model["daily"]
globe_daily = globe_daily[...,:360]

indices = np.load(tracer_path)
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
bac = tracer_indices["bac"][0]
phy = tracer_indices["phyto"][0]
zoo = tracer_indices["zoo"][0]
don = tracer_indices["don"][0]
pon = tracer_indices["pon"][0]


xlabel = ['J','A','J','O','J']
xticks = [0,90,180,270,360]
days = np.linspace(0,359,360)

fig, axs = plt.subplots(2,4,figsize=(20,10), sharex=True)

axs[0,0].plot(days,globe_daily[no3,0],'-b',label='GLOBE')
axs[0,0].plot(days,fasham_daily[no3],'-.k',label='Fasham (1990) - FABM')
axs[0,0].set_title("(a) Nitrate")
axs[0,0].set_xlim([0,360])
axs[0,0].set_xlabel("Days")
axs[0,0].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[0,1].plot(days,globe_daily[bac,0],'-b',label='GLOBE')
axs[0,1].plot(days,fasham_daily[bac],'-.k',label='Fasham (1990) - FABM')
axs[0,1].set_title("(b) Bacteria")
axs[0,1].set_xlim([0,360])
axs[0,1].set_xlabel("Days")
axs[0,1].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[0,2].plot(days,globe_daily[phy,0],'-b',label='GLOBE')
axs[0,2].plot(days,fasham_daily[phy],'-.k',label='Fasham (1990) - FABM')
axs[0,2].set_title("(c) Phytoplankton")
axs[0,2].set_xlim([0,360])
axs[0,2].set_xlabel("Days")
axs[0,2].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[0,3].plot(days,globe_daily[zoo,0],'-b',label='GLOBE')
axs[0,3].plot(days,fasham_daily[zoo],'-.k',label='Fasham (1990) - FABM')
axs[0,3].set_title("(d) Zooplankton")
axs[0,3].set_xlim([0,360])
axs[0,3].set_xlabel("Days")
axs[0,3].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[1,0].plot(days,globe_daily[nh4,0],'-b',label='GLOBE')
axs[1,0].plot(days,fasham_daily[nh4],'-.k',label='Fasham (1990) - FABM')
axs[1,0].set_title("(e) Ammonium")
axs[1,0].set_xlim([0,360])
axs[1,0].set_xlabel("Days")
axs[1,0].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[1,1].plot(days,globe_daily[don,0],'-b',label='GLOBE')
axs[1,1].plot(days,fasham_daily[don],'-.k',label='Fasham (1990) - FABM')
axs[1,1].set_title("(f) Dissolved Organic Nitrogen")
axs[1,1].set_xlim([0,360])
axs[1,1].set_xlabel("Days")
axs[1,1].set_ylabel("mmol $N$ ${m^{-3}}$")

axs[1,2].plot(days,globe_daily[pon,0],'-b',label='GLOBE')
axs[1,2].plot(days,fasham_daily[pon],'-.k',label='Fasham (1990) - FABM')
axs[1,2].set_title("(g) Detritus")
axs[1,2].set_xlim([0,360])
axs[1,2].set_xlabel("Days")
axs[1,2].set_ylabel("mmol $N$ ${m^{-3}}$")

handles,labels = axs[0,0].get_legend_handles_labels()
axs[1,3].axis('off')
axs[1,3].legend(handles,labels,loc='center',frameon=True)

folder = os.getcwd() + '/tests/fabm'
fig.tight_layout()
save_loc = os.path.join(folder + '/figures','fasham0d_seasonal.jpg')
plt.savefig(save_loc)


