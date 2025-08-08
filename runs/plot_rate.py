# Code to plot modulation simulation results as rate vs year

# Import modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import sys

# Check specie
if len(sys.argv) < 2:
   print("Error: Specie must be specified.")
   print("Accepted species: 'helium', 'hydrogen', or 'electrons'.")
   exit(1)
elif sys.argv[1] == "helium":
   specie_label = "He"
elif sys.argv[1] == "hydrogen":
   specie_label = "H"
elif sys.argv[1] == "electrons":
   specie_label = "e"
else:
   print("Error: Unrecognized specie.")
   print("Accepted species: 'helium', 'hydrogen', or 'electrons'.")
   exit(1)
   
print("Plotting results for {:s}.".format(sys.argv[1]))

# Import simulation data
file_names = [
              "../results/HS_mod_spec_{:s}/HS_mod_parker_integ_spec.dat".format(specie_label),
              # "../results/HS_mod_spec_{:s}_no_MAX/HS_mod_parker_integ_spec.dat".format(specie_label),
              # "../results/HS_mod_spec_{:s}_no_UHS/HS_mod_parker_integ_spec.dat".format(specie_label),
              # "../results/HS_mod_spec_{:s}_no_MAX_or_UHS/HS_mod_parker_integ_spec.dat".format(specie_label),
              # "../results/HS_mod_spec_{:s}_no_drift/HS_mod_parker_integ_spec.dat".format(specie_label),
              ]
labels = [
          "full",
          # "no MAX",
          # "no UHS",
          # "no UHS or MAX",
          # "no drift",
          ]
markers = ["o","^","s","X","D","P"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
          "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]
num_data_files = len(file_names)
data = [None for _ in range(num_data_files)]
for file in range(num_data_files):
   data[file] = np.loadtxt(file_names[file])

# Import Voyager trajectory and make interpolation function to convert radius to year
V2_path_file_name = "data/V2_RTP_HGI.dat"
V2_traj = np.loadtxt(V2_path_file_name, skiprows=1)
V2_year = V2_traj[:,0] + (V2_traj[:,1] - 1) / (365 + (4 - V2_traj[:,0] % 4) // 4)
V2_path = V2_traj[:,2]
rad2year = interp1d(V2_path, V2_year)

# Import Voyager and HCS latitude along trajectory
V2_lat_file_name = "../results/V2_vs_HCS_lat.dat"
V2_lat = np.loadtxt(V2_lat_file_name)
UHS_seg = []
HS_type_crossings = []
for i in range(np.size(V2_lat,0)-1):
   if (V2_lat[i,4]*V2_lat[i+1,4] < 0):
      HS_type_crossings.append(0.5 * (V2_lat[i,1]+V2_lat[i+1,1]))
if V2_lat[0,4] == -1:
   HS_type_crossings.insert(0, V2_lat[0,1])
if V2_lat[-1,4] == -1:
   HS_type_crossings.append(V2_lat[-1,1])
for i in range(len(HS_type_crossings)//2):
   UHS_seg.append([HS_type_crossings[2*i], HS_type_crossings[2*i+1]])

# Import Voyager data
V2_rate_file_name = "data/V2_{:s}_rate.dat".format(specie_label)
V2_rate = np.loadtxt(V2_rate_file_name)

# Separate data based on background threshold
threshold = 0.5
cut_idx = np.size(V2_rate[:,0])
for pt in range(np.size(V2_rate[:,0])-1,-1,-1):
   if threshold * V2_rate[pt,1] < V2_rate[pt,2]:
      cut_idx = pt+1
      break

# Plot data
fig = plt.figure(figsize=(12, 8), layout='tight')
ax1 = fig.add_subplot(111, projection='rectilinear')
ax2 = ax1.twiny()

if sys.argv[1] == "electrons":
# Correct signal by subtracting background
   V2_rate[:,1] = V2_rate[:,1] - V2_rate[:,2]
   ax1.semilogy(V2_rate[:cut_idx+1,0], V2_rate[:cut_idx+1,1], color="c", linewidth=2, zorder=0, label="Observations (bkg > {:.0f}%)".format(threshold*100))
   if cut_idx < np.size(V2_rate[:,0]):
      ax1.semilogy(V2_rate[cut_idx:,0], V2_rate[cut_idx:,1], color="k", linewidth=2, zorder=0, label="Observations (bkg < {:.0f}%)".format(threshold*100))
else:
   ax1.plot(V2_rate[:cut_idx+1,0], V2_rate[:cut_idx+1,1], color="c", linewidth=2, zorder=0, label="Observations (bkg > {:.0f}%)".format(threshold*100))
   if cut_idx < np.size(V2_rate[:,0]):
      ax1.plot(V2_rate[cut_idx:,0], V2_rate[cut_idx:,1], color="k", linewidth=2, zorder=0, label="Observations (bkg < {:.0f}%)".format(threshold*100))
for seg in range(len(UHS_seg)):
   ax1.axvspan(UHS_seg[seg][0], UHS_seg[seg][1], alpha=0.25, color='red')
for file in range(num_data_files):
   ax1.scatter(rad2year(data[file][:,0]), data[file][:,1], s=80, marker=markers[file], c=colors[file], label=labels[file])

# Bottom axis
ax1.set_xlabel('Year', fontsize=20)
ax1.set_ylabel(specie_label + " Rate (s$^{-1}$)", fontsize=20)
for i in range(np.size(V2_path)):
   if V2_year[i] < 2007.5:
      idx_left = i
   else:
      break
for i in range(idx_left, np.size(V2_path)):
   if V2_year[i] < 2018.83:
      idx_right = i
   else:
      break
ax1.set_xlim(V2_year[idx_left], V2_year[idx_right])
ax1.tick_params(labelsize=20)
ax1.legend(loc=4, fontsize=20)

# Top axis
ax2.set_xlim(V2_path[idx_left], V2_path[idx_right])
ax2.set_xticks([84.0, 92.0, 100.0, 108.0, 116.0])
ax2.tick_params(labelsize=20)
ax2.set_xlabel("Radial Distance (au)", fontsize=20)

# Vertical lines
y_bot, y_top = ax1.get_ylim()
a_bot = 0.1
a_top = 1.0 - a_bot
if sys.argv[1] == "electrons":
   y_mid = np.exp(a_bot * np.log(y_bot) + a_top * np.log(y_top))
else:
   y_mid = a_bot * y_bot + a_top * y_top
ax1.annotate("TS", (2007.8,y_mid), fontsize=24)
ax1.axvline(2007.67, color='k', linestyle='--', linewidth=2)
ax1.annotate("MAX", (2012.05,y_mid), color='m', fontsize=24)
ax1.axvline(2013.05, color='m', linestyle=':', linewidth=2)

plt.show()
plt.close(fig)
