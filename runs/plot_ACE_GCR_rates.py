# Code to plot ACE CRIS rates at 1au

import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
import time

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

# Import data
data = np.loadtxt("data/ACE_GCR_rate.dat", skiprows=1)
year_float = np.zeros(np.size(data,0))

# Convert (year, doy) to year-float
for i in range(np.size(data,0)):
   date0 = dt(year=int(data[i,0]),month=1,day=1,hour=0,minute=0,second=0)
   date = date0 + td(days=int(data[i,1])-1,hours=int(data[i,2]),minutes=0)
   year_float[i] = toYearFraction(date)

rate = data[:,3] # get hourly rate
Imputate(year_float, rate, [0.000]) # Interpolate 0.0 values
rate_smooth = savgol_filter(rate, 24*52, 0) # 52 day average
rate_smoother = savgol_filter(rate, 24*365, 0) # 1 year average
rate_detrend = rate_smooth - rate_smoother # find detrended rate

# Plot
fig = plt.figure(figsize=(12, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.plot(year_float, rate_detrend, label="1 au", linewidth=2)
ax.set_xlabel('Year', fontsize=20)
ax.set_ylabel("CRIS E9 rate (s$^{-1}$)", fontsize=20)
ax.set_xlim(2005,2015)
ax.tick_params(labelsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.legend(fontsize=16)

plt.show()
plt.close(fig)