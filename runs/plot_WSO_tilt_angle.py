# Code to fit WSO observations

# Import modules
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime as dt
import time

# Cubic stretch function
def CubicStretch(t, d):
   return t * ((1.0 - d) * (t / np.pi)**2 + 3.0 * (d - 1.0) * (t / np.pi) + 3.0 - 2.0 * d)

# Periodic function
def HCS(t, p0, A, T, t0, d):
   arg = 4.0 * np.pi * (t - t0) / T
   arg_mod = np.floor(arg / (2.0 * np.pi))
   arg_rem = CubicStretch(arg - arg_mod * 2.0 * np.pi, d)
   return p0 + A * np.cos(arg_rem)

# Date to year (float) function
def toYearFraction(date):
   def sinceEpoch(date): # returns seconds since epoch
      return time.mktime(date.timetuple())
   s = sinceEpoch

   year = date.year
   startOfThisYear = dt(year=year, month=1, day=1)
   startOfNextYear = dt(year=year+1, month=1, day=1)

   yearElapsed = s(date) - s(startOfThisYear)
   yearDuration = s(startOfNextYear) - s(startOfThisYear)
   fraction = yearElapsed/yearDuration

   return date.year + fraction

# Import data
CRi = 1642
CRf = 2285
year_float = []
CR = []
Rav = []
Rn = []
Rs = []
Lav = []
Ln = []
Ls = []
LRav = []
LRn = []
LRs = []
aL = 0.2
aR = 1.0 - aL

file = open("data/WSO_tilt_angle.dat", "r")
# Header lines
for row in range(1):
   line = file.readline()
# Data
line = file.readline()
while line:
   data_str = line.split()
   CRn = int(data_str[1]) # Carrington rotation number
   if (CRi <= CRn <= CRf):
      CR.append(CRn)
      date = dt.strptime(data_str[2],"%Y:%m:%d")
      year_float.append(toYearFraction(date))
      Rav.append(float(data_str[4]))
      Rn.append(float(data_str[5]))
      Rs.append(-float(data_str[6]))
      Lav.append(float(data_str[7]))
      Ln.append(float(data_str[8]))
      Ls.append(-float(data_str[9]))
      LRav.append(aL * Lav[-1] + aR * Rav[-1])
      LRn.append(aL * Ln[-1] + aR * Rn[-1])
      LRs.append(aL * Ls[-1] + aR * Rs[-1])
   line = file.readline()
file.close()
year_float = np.array(year_float)

# Best fit periodic function
opt_params, cov = curve_fit(HCS, year_float, Rav,
                            p0=[40.0, 35.0, 22.0, 2002.0, 0.5], maxfev=10000,
                            bounds=([40.0, 35.0, 21.9, 2002.0, 0.0], [40.1, 35.1, 22.0, 2002.1, 1.0]))
print("Optimal Fit Parameters:", opt_params)

# Slice
WSO_slice = slice(0,-1,1)

# Plot
fig = plt.figure(figsize=(12, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.plot(year_float, LRs, linestyle="", marker="o", label="WSO south", markersize=8)
ax.plot(year_float[WSO_slice], LRs[WSO_slice], label="slice", linewidth=4)
ax.set_xlabel('Year', fontsize=30)
ax.set_ylabel('Computed HCS Tilt Angle ($^\\circ$)', fontsize=30)
ax.set_ylim(0.0,90.0)
ax.tick_params(labelsize=30)
ax.legend(fontsize=20)

plt.show()
plt.close(fig)

# Output slices
indices = np.arange(CRf - CRi + 1)

file = open("data/WSO_tilt_angle_slice_LRs.dat", "w")
file.write("{:d}\n".format(len(indices[WSO_slice])))
for i in indices[WSO_slice]:
   file.write("{:12.6f}{:12.6f}\n".format(year_float[i], LRs[i]))
file.close()
