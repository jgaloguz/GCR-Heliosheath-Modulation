# Code to plot V2 CRS rates in HS

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.signal import savgol_filter
from datetime import datetime as dt
import time

# Import data
H_rate = np.loadtxt("data/V2_H_rate.dat")
He_rate = np.loadtxt("data/V2_He_rate.dat")
e_rate = np.loadtxt("data/V2_e_rate.dat")

# Plot
fig = plt.figure(figsize=(12, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

axcolor = 'tab:red'
ax.plot(H_rate[:,0], H_rate[:,1]-H_rate[:,2], color=axcolor, label="H", linewidth=2)
ax.plot(He_rate[:,0], 10.0*(He_rate[:,1]-He_rate[:,2]), color=axcolor, label="He (x10)", linestyle=":", linewidth=2)
ax.set_xlabel('Year', fontsize=20)
ax.set_ylabel("Ion rate (s$^{-1}$)", fontsize=20, color=axcolor)
ax.axvline(2007.67, color='k', linestyle='--', linewidth=2)
ax.annotate("TS", (2007.8,0.5), fontsize=24)
ax.axvline(2018.83, color='k', linestyle='--', linewidth=2)
ax.annotate("HP", (2019.0,0.5), fontsize=24)
ax.set_xlim(2007,2020)
ax.tick_params(labelsize=20)
ax.tick_params(axis='x', which='minor', labelsize=20)
ax.tick_params(axis='y', labelcolor=axcolor)
ax.legend(fontsize=16)

axx = ax.twinx()
axxcolor = 'tab:blue'
axx.semilogy(e_rate[:,0], e_rate[:,1]-e_rate[:,2], color=axxcolor, label="e", linewidth=2)
axx.set_ylabel("Electron rate (s$^{-1}$)", fontsize=20, color=axxcolor)
axx.tick_params(labelsize=20)
axx.tick_params(axis='y', labelcolor=axxcolor)

plt.show()
plt.close(fig)
