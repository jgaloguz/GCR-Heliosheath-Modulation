# Code to plot flow and field along V2 trajectory

# Import libraries
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
import time
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# Function to convert date (year, doy, minute) to year (float)
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

# Function to fill in missing values of data by interpolating bounding valid data
def Imputate(time_array, data_array, missing_value):
   N = np.size(data_array)
# Fix first value if necessary by setting it equal to the first non-zero value
   i = 0
   while data_array[i] == missing_value:
      i = i + 1
   data_array[0] = data_array[i]
# Fix last value if necessary by setting it equal to the last non-zero value
   i = N-1
   while data_array[i] == missing_value:
      i = i - 1
   data_array[N-1] = data_array[i]
# Interpolate missing_values between first and last
   i = 1
   while i < N-1:
# If value is missing
      if data_array[i] == missing_value:
# Get last value not missing (just previous value)
         k1 = i-1
         j = i+1
         while data_array[j] == missing_value:
            j = j + 1
# Get next value not missing (guaranteed to exist)
         k2 = j
# Interpolate
         for k in range(k1+1,k2):
            data_array[k] = (data_array[k1] * (time_array[k2] - time_array[k]) + data_array[k2] * (time_array[k] - time_array[k1])) / (time_array[k2] - time_array[k1])
         i = j
      i = i + 1

# Import data with gaps
bad_value = [999.999, 9999.9, 99.99999] # Bad value flags for each input
labels = ["nT", "km/s", "cm$^{-3}$"] # label for each input
data = np.loadtxt("data/V2_daily_flow_speed_mag_field.dat")
year = data[:,0]
day = data[:,1]
V2_r = data[:,3]
Bmag = data[:,4]
Vmag = data[:,5]
Np = data[:,6]
n_days = np.size(year)
year_float = np.zeros(n_days)
for t in range(n_days):
   date0 = dt(year=int(year[t]), month=1, day=1)
   date = date0 + td(days=int(day[t]))
   year_float[t] = toYearFraction(date)

# Fill in bad values with linear interpolation
Imputate(year_float, Bmag, bad_value[0])
Imputate(year_float, Vmag, bad_value[1])
Imputate(year_float, Np, bad_value[2])

# Smooth arrays
Bmag = savgol_filter(Bmag, 30, 0)
Vmag = savgol_filter(Vmag, 30, 0)
Np = savgol_filter(Np, 30, 0)

# Import simulation data
data_sim = np.loadtxt("../results/V2_flow_speed_mag_field_sim.dat")
V2_r_sim = data_sim[:,0] # au
year_float_sim = data_sim[:,1] # yr
Vmag_sim = data_sim[:,2] * 1.0e-5 # cm/s -> km/s
Bmag_sim = data_sim[:,3] * 1.0e5  # G    -> nT

# Plot
fig = plt.figure(figsize=(14, 10), layout='tight')
ax1 = fig.add_subplot(211, projection='rectilinear')

ax1.plot(year_float, Vmag, 'b-', linewidth=2, label='data')
ax1.plot(year_float_sim, Vmag_sim, 'r-', linewidth=2, label='sim')
ax1.set_xlabel('Year', fontsize=20)
ax1.set_ylabel('$|V|$ (km/s)', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.set_xlim(2007,2018.84)
ax1.set_ylim(0.0,410.0)
ax1.axvline(2007.65, color='k', linestyle='--', linewidth=2)
ax1.legend(fontsize=20)

ax2 = fig.add_subplot(212, projection='rectilinear')

ax2.plot(year_float, Bmag, 'b-', linewidth=2, label='data')
ax2.plot(year_float_sim, Bmag_sim, 'r-', linewidth=2, label='sim')
ax2.set_xlabel('Year', fontsize=20)
ax2.set_ylabel('$|B|$ (nT)', fontsize=20)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlim(2007,2018.83)
ax2.set_ylim(0.0,0.5)
ax2.annotate("TS", (2007.8,0.35), fontsize=24)
ax2.axvline(2007.65, color='k', linestyle='--', linewidth=2)
ax2.annotate("MAX", (2012.25,0.35), color='m', fontsize=24)
ax2.axvline(2013.05, color='m', linestyle=':', linewidth=2)

plt.show()
plt.close(fig)