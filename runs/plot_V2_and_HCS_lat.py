# Code to plot V2 latitude vs advected HCS latitude

# Import modules
import numpy as np
import matplotlib.pyplot as plt

# Load data
lat_data = np.loadtxt("../results/V2_vs_HCS_lat.dat")
UHS_seg = []
HS_type_crossings = []
for i in range(np.size(lat_data,0)-1):
   if (lat_data[i,4]*lat_data[i+1,4] < 0):
      HS_type_crossings.append(0.5 * (lat_data[i,1]+lat_data[i+1,1]))
if lat_data[0,4] == -1:
   HS_type_crossings.insert(0, lat_data[0,1])
if lat_data[-1,4] == -1:
   HS_type_crossings.append(lat_data[-1,1])
for i in range(len(HS_type_crossings)//2):
   UHS_seg.append([HS_type_crossings[2*i], HS_type_crossings[2*i+1]])

# Plot
fig = plt.figure(figsize=(12, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.plot(lat_data[:,2], lat_data[:,1], label="V2", linewidth=3)
ax.plot(lat_data[:,3], lat_data[:,1], label="WSO", linewidth=3)
for seg in range(len(UHS_seg)):
   ax.axhspan(UHS_seg[seg][0], UHS_seg[seg][1], alpha=0.25, color='red')
ax.set_xlabel('Southern Latitude ($^\\circ$)', fontsize=30)
ax.set_ylabel('Year', fontsize=30)
ax.set_xlim(0.0,80.0)
ax.set_ylim(lat_data[0,1],lat_data[-1,1])
ax.tick_params(labelsize=30)
ax.legend(loc=2, fontsize=20)

plt.show()
plt.close(fig)