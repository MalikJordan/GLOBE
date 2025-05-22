import numpy as np
from matplotlib import pyplot as plt
import os

# Folder path
folder = os.getcwd() + '/tests/bfm17/0522'

# Load solution
solution = np.load("bfm17-0522-fdm.npz", allow_pickle=True)
# conc = solution["conc"]     # concentratrion matrix
conc = solution["concentration"]     # concentratrion matrix
time = solution["time"]     # time array

# Load tracer indices
indices = np.load("tracer_indices_bfm17-0522.npz")
tracer_indices = {}
for file in indices.files:
    tracer_indices[file] = list(indices[file])

o2 = tracer_indices["o2"][0]
no3 = tracer_indices["no3"][0]
nh4 = tracer_indices["nh4"][0]
po4 = tracer_indices["po4"][0]
pc = tracer_indices["phyto1"][0]
pn = tracer_indices["phyto1"][1]
pp = tracer_indices["phyto1"][2]
pl = tracer_indices["phyto1"][3]
zc = tracer_indices["zoo1"][0]
zn = tracer_indices["zoo1"][1]
zp = tracer_indices["zoo1"][2]
domc = tracer_indices["dom1"][0]
domn = tracer_indices["dom1"][1]
domp = tracer_indices["dom1"][2]
pomc = tracer_indices["pom1"][0]
pomn = tracer_indices["pom1"][1]
pomp = tracer_indices["pom1"][2]

# Daily averages for concentrations
iters = np.linspace(0,len(time)-1,len(time))
iters_per_day = 86400/360       # change this depending on timestep
conc_avg = np.zeros((conc.shape[0], int(np.floor((len(iters))/iters_per_day))))

day = 0
for i in range(0,len(iters)):
    conc_avg[:,day] += conc[:,i]
    if i%int(iters_per_day) == 0 and i != 0:
        conc_avg[:,day] = conc_avg[:,day]/iters_per_day

        day += 1
        if day == 730:
            break

# Convert time array from seconds to months
sec_mon = 60 * 60 * 24 * 30
time = time/sec_mon
days = np.linspace(0,729,730)

# Ticks and labels for plots
# xticks = [0,2,4,6,8,10,12,14,16,18,20,22,24]
xlabel = ['J','M','M','J','S','N','J','M','M','J','S','N','J']

xticks = [0,60,120,180,240,300,360,420,480,540,600,660,720]

# Create plots
# Nutrients
fig, axs = plt.subplots(2,2,figsize=(12,12))

# axs[0,0].plot(time,conc[o2])
axs[0,0].plot(days,conc_avg[o2])
axs[0,0].set_title("Oxygen")
axs[0,0].set_xlabel("Time [months]")
axs[0,0].set_xticks(xticks,xlabel)
# axs[0,0].set_xlim([0,24])
axs[0,0].set_xlim([0,730])
axs[0,0].set_ylabel("mmol O $\mathregular{m^{-3}}$")

# axs[0,1].plot(time,conc[no3])
axs[0,1].plot(days,conc_avg[no3])
axs[0,1].set_title("Nitrate")
axs[0,1].set_xlabel("Time [months]")
axs[0,1].set_xticks(xticks,xlabel)
# axs[0,1].set_xlim([0,24])
axs[0,1].set_xlim([0,730])
axs[0,1].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1,0].plot(time,conc[nh4])
axs[1,0].plot(days,conc_avg[nh4])
axs[1,0].set_title("Ammonium")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
# axs[1,0].set_xlim([0,24])
axs[1,0].set_xlim([0,730])
axs[1,0].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1,1].plot(time,conc[po4])
axs[1,1].plot(days,conc_avg[po4])
axs[1,1].set_title("Phosphate")
axs[1,1].set_xlabel("Time [months]")
axs[1,1].set_xticks(xticks,xlabel)
# axs[1,1].set_xlim([0,24])
axs[1,1].set_xlim([0,730])
axs[1,1].set_ylabel("mmol P $\mathregular{m^{-3}}$")

fig.suptitle("Nutrients")
fig.tight_layout()
nut = os.path.join(folder,"nutrients.jpg")
plt.savefig(nut)

# Phytoplankton
fig, axs = plt.subplots(2,2,figsize=(12,12))

# axs[0,0].plot(time,conc[pc])
axs[0,0].plot(days,conc_avg[pc])
axs[0,0].set_title("Carbon")
axs[0,0].set_xlabel("Time [months]")
axs[0,0].set_xticks(xticks,xlabel)
# axs[0,0].set_xlim([0,24])
axs[0,0].set_xlim([0,730])
axs[0,0].set_ylabel("mg C $\mathregular{m^{-3}}$")

# axs[0,1].plot(time,conc[pn])
axs[0,1].plot(days,conc_avg[pn])
axs[0,1].set_title("Nitrogen")
axs[0,1].set_xlabel("Time [months]")
axs[0,1].set_xticks(xticks,xlabel)
# axs[0,1].set_xlim([0,24])
axs[0,1].set_xlim([0,730])
axs[0,1].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[1,0].plot(time,conc[pp])
axs[1,0].plot(days,conc_avg[pp])
axs[1,0].set_title("Phosphorus")
axs[1,0].set_xlabel("Time [months]")
axs[1,0].set_xticks(xticks,xlabel)
# axs[1,0].set_xlim([0,24])
axs[1,0].set_xlim([0,730])
axs[1,0].set_ylabel("mmol P $\mathregular{m^{-3}}$")

# axs[1,1].plot(time,conc[pl])
axs[1,1].plot(days,conc_avg[pl])
axs[1,1].set_title("Chlorophyll-a")
axs[1,1].set_xlabel("Time [months]")
axs[1,1].set_xticks(xticks,xlabel)
# axs[1,1].set_xlim([0,24])
axs[1,1].set_xlim([0,730])
axs[1,1].set_ylabel("mg Chl-a $\mathregular{m^{-3}}$")

fig.suptitle("Phytoplankton")
fig.tight_layout()
phyto = os.path.join(folder,"phytoplankton.jpg")
plt.savefig(phyto)

# Zooplankton
fig, axs = plt.subplots(1,3,figsize=(15,5))

# axs[0].plot(time,conc[zc])
axs[0].plot(days,conc_avg[zc])
axs[0].set_title("Carbon")
axs[0].set_xlabel("Time [months]")
axs[0].set_xticks(xticks,xlabel)
# axs[0].set_xlim([0,24])
axs[0].set_xlim([0,730])
axs[0].set_ylabel("mg C $\mathregular{m^{-3}}$")

# axs[0].plot(time,conc[zn])
axs[1].plot(days,conc_avg[zn])
axs[1].set_title("Nitrogen")
axs[1].set_xlabel("Time [months]")
axs[1].set_xticks(xticks,xlabel)
# axs[1].set_xlim([0,24])
axs[1].set_xlim([0,730])
axs[1].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[2].plot(time,conc[zp])
axs[2].plot(days,conc_avg[zp])
axs[2].set_title("Phosphorus")
axs[2].set_xlabel("Time [months]")
axs[2].set_xticks(xticks,xlabel)
# axs[2].set_xlim([0,24])
axs[2].set_xlim([0,730])
axs[2].set_ylabel("mmol P $\mathregular{m^{-3}}$")

fig.suptitle("Zooplankton")
fig.tight_layout()
zoo = os.path.join(folder,"zooplankton.jpg")
plt.savefig(zoo)

# DOM
fig, axs = plt.subplots(1,3,figsize=(15,5))

# axs[0].plot(time,conc[domc])
axs[0].plot(days,conc_avg[domc])
axs[0].set_title("Carbon")
axs[0].set_xlabel("Time [months]")
axs[0].set_xticks(xticks,xlabel)
# axs[0].set_xlim([0,24])
axs[0].set_xlim([0,730])
axs[0].set_ylabel("mg C $\mathregular{m^{-3}}$")

# axs[1].plot(time,conc[domn])
axs[1].plot(days,conc_avg[domn])
axs[1].set_title("Nitrogen")
axs[1].set_xlabel("Time [months]")
axs[1].set_xticks(xticks,xlabel)
# axs[1].set_xlim([0,24])
axs[1].set_xlim([0,730])
axs[1].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[2].plot(time,conc[domp])
axs[2].plot(days,conc_avg[domp])
axs[2].set_title("Phosphorus")
axs[2].set_xlabel("Time [months]")
axs[2].set_xticks(xticks,xlabel)
# axs[2].set_xlim([0,24])
axs[2].set_xlim([0,730])
axs[2].set_ylabel("mmol P $\mathregular{m^{-3}}$")

fig.suptitle("Dissolved Organic Matter")
fig.tight_layout()
dom = os.path.join(folder,"dom.jpg")
plt.savefig(dom)

# DOM
fig, axs = plt.subplots(1,3,figsize=(15,5))

# axs[0].plot(time,conc[pomc])
axs[0].plot(days,conc_avg[pomc])
axs[0].set_title("Carbon")
axs[0].set_xlabel("Time [months]")
axs[0].set_xticks(xticks,xlabel)
# axs[0].set_xlim([0,24])
axs[0].set_xlim([0,730])
axs[0].set_ylabel("mg C $\mathregular{m^{-3}}$")

# axs[1].plot(time,conc[pomn])
axs[1].plot(days,conc_avg[pomn])
axs[1].set_title("Nitrogen")
axs[1].set_xlabel("Time [months]")
axs[1].set_xticks(xticks,xlabel)
# axs[1].set_xlim([0,24])
axs[1].set_xlim([0,730])
axs[1].set_ylabel("mmol N $\mathregular{m^{-3}}$")

# axs[2].plot(time,conc[pomp])
axs[2].plot(days,conc_avg[pomp])
axs[2].set_title("Phosphorus")
axs[2].set_xlabel("Time [months]")
axs[2].set_xticks(xticks,xlabel)
# axs[2].set_xlim([0,24])
axs[2].set_xlim([0,730])
axs[2].set_ylabel("mmol P $\mathregular{m^{-3}}$")

fig.suptitle("Particulate Organic Matter")
fig.tight_layout()
pom = os.path.join(folder,"pom.jpg")
plt.savefig(pom)