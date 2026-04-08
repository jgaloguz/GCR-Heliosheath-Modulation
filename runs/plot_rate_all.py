# Code to plot modulation simulation results as rate vs year

# Import modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import sys

# No terminal inputs

# Import simulation data
panel = "left" # Plot left or right panel of main simulation results figure
if panel == "left":
   markers = ["o","s","^","X"]
   colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
   file_names_He = [
                 "../results/HS_mod_spec_He_R1/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R2/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R3/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R4/HS_mod_parker_integ_spec.dat",
                 ]
   file_names_H = [
                 "../results/HS_mod_spec_H_R1/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R2/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R3/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R4/HS_mod_parker_integ_spec.dat",
                 ]
   file_names_e = [
                 "../results/HS_mod_spec_e_R1/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R2/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R3/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R4/HS_mod_parker_integ_spec.dat",
                 ]
   labels = ["R1", "R2", "R3", "R4"]
else:
   file_names_He = [
                 "../results/HS_mod_spec_He_R5/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R6/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R7/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_He_R8/HS_mod_parker_integ_spec.dat",
                 ]
   file_names_H = [
                 "../results/HS_mod_spec_H_R5/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R6/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R7/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_H_R8/HS_mod_parker_integ_spec.dat",
                 ]
   file_names_e = [
                 "../results/HS_mod_spec_e_R5/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R6/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R7/HS_mod_parker_integ_spec.dat",
                 "../results/HS_mod_spec_e_R8/HS_mod_parker_integ_spec.dat",
                 ]
   labels = ["R5", "R6", "R7", "R8"]
   markers = ["D","P","H","*"]
   colors = ["tab:purple", "tab:brown", "tab:cyan", "tab:olive"]
num_data_files_He = len(file_names_He)
data_He = [None for _ in range(num_data_files_He)]
for file in range(num_data_files_He):
   data_He[file] = np.loadtxt(file_names_He[file])
num_data_files_H = len(file_names_H)
data_H = [None for _ in range(num_data_files_H)]
for file in range(num_data_files_H):
   data_H[file] = np.loadtxt(file_names_H[file])
num_data_files_e = len(file_names_e)
data_e = [None for _ in range(num_data_files_e)]
for file in range(num_data_files_e):
   data_e[file] = np.loadtxt(file_names_e[file])

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
V2_rate_file_name_He = "data/V2_He_rate.dat"
V2_rate_file_name_H = "data/V2_H_rate.dat"
V2_rate_file_name_e = "data/V2_e_rate.dat"
V2_rate_He = np.loadtxt(V2_rate_file_name_He)
V2_rate_H = np.loadtxt(V2_rate_file_name_H)
V2_rate_e = np.loadtxt(V2_rate_file_name_e)

# Separate data based on background threshold
n_V2 = np.size(V2_rate_He, 0)
threshold = 0.5
cut_idx_He = [0, n_V2-1]
flag = True
for pt in range(n_V2-1,-1,-1):
   if flag and threshold * V2_rate_He[pt,1] < V2_rate_He[pt,2]:
      cut_idx_He.insert(1, pt+1)
      flag = False
   elif not flag and threshold * V2_rate_He[pt,1] > V2_rate_He[pt,2]:
      cut_idx_He.insert(1, pt+1)
      flag = True
n_V2 = np.size(V2_rate_H, 0)
threshold = 0.5
cut_idx_H = [0, n_V2-1]
flag = True
for pt in range(n_V2-1,-1,-1):
   if flag and threshold * V2_rate_H[pt,1] < V2_rate_H[pt,2]:
      cut_idx_H.insert(1, pt+1)
      flag = False
   elif not flag and threshold * V2_rate_H[pt,1] > V2_rate_H[pt,2]:
      cut_idx_H.insert(1, pt+1)
      flag = True
n_V2 = np.size(V2_rate_e, 0)
threshold = 0.5
cut_idx_e = [0, n_V2-1]
flag = True
for pt in range(n_V2-1,-1,-1):
   if flag and threshold * V2_rate_e[pt,1] < V2_rate_e[pt,2]:
      cut_idx_e.insert(1, pt+1)
      flag = False
   elif not flag and threshold * V2_rate_e[pt,1] > V2_rate_e[pt,2]:
      cut_idx_e.insert(1, pt+1)
      flag = True

# Correct electron signal by subtracting background
V2_rate_e[:,1] = V2_rate_e[:,1] - V2_rate_e[:,2]

# Plot data
a_bot = 0.1
a_top = 1.0 - a_bot
fig = plt.figure(figsize=(12, 24), layout='tight')

ax1 = fig.add_subplot(311, projection='rectilinear')
ax2 = ax1.twiny()

# Iterate over observation segments based on background threshold
seg_color = "k"
for seg in range(len(cut_idx_He)-1, 0, -1):
   obs, = ax1.plot(V2_rate_He[cut_idx_He[seg-1]:cut_idx_He[seg]+1,0], V2_rate_He[cut_idx_He[seg-1]:cut_idx_He[seg]+1,1],
                    color=seg_color, linewidth=2, zorder=0)
   if seg_color == "k":
      seg_color = "c"
   else:
      seg_color = "k"
   if seg == len(cut_idx_He)-1:
      obs.set_label("V2 (bkg < {:.0f}%)".format(threshold*100))
   elif seg == len(cut_idx_He)-2:
      obs.set_label("V2 (bkg > {:.0f}%)".format(threshold*100))

# Color background of image depending on UHS or SHS
for seg in range(len(UHS_seg)):
   ax1.axvspan(UHS_seg[seg][0], UHS_seg[seg][1], alpha=0.25, color='red')

# Plot simulation results
for file in range(num_data_files_He):
   ax1.scatter(rad2year(data_He[file][:,0]), data_He[file][:,1], s=80, marker=markers[file], c=colors[file], label=labels[file])

# Bottom axis
ax1.set_xlabel('Year', fontsize=20)
ax1.set_ylabel("He$^{2\\!+}$ Rate (s$^{-1}$)", fontsize=20)
for i in range(np.size(V2_path)):
   if V2_year[i] < 2007.0:
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
ax1.set_ylim(0.02,0.045)
ax1.legend(fontsize=20)

# Top axis
ax2.set_xlim(V2_path[idx_left], V2_path[idx_right])
ax2.set_xticks([84.0, 92.0, 100.0, 108.0, 116.0])
ax2.tick_params(labelsize=20)
ax2.set_xlabel("Radial Distance (au)", fontsize=20)

# Vertical lines
y_bot, y_top = ax1.get_ylim()
y_mid = a_bot * y_bot + a_top * y_top
ax1.annotate("TS", (2007.8,y_mid), fontsize=24)
ax1.axvline(2007.67, color='k', linestyle='--', linewidth=2)
ax1.annotate("MAX", (2012.05,y_mid), color='m', fontsize=24)
ax1.axvline(2013.05, color='m', linestyle=':', linewidth=2)

ax1 = fig.add_subplot(312, projection='rectilinear')
ax2 = ax1.twiny()

# Iterate over observation segments based on background threshold
seg_color = "k"
for seg in range(len(cut_idx_H)-1, 0, -1):
   obs, = ax1.plot(V2_rate_H[cut_idx_H[seg-1]:cut_idx_H[seg]+1,0], V2_rate_H[cut_idx_H[seg-1]:cut_idx_H[seg]+1,1],
                    color=seg_color, linewidth=2, zorder=0)
   if seg_color == "k":
      seg_color = "c"
   else:
      seg_color = "k"
   if seg == len(cut_idx_H)-1:
      obs.set_label("V2 (bkg < {:.0f}%)".format(threshold*100))
   elif seg == len(cut_idx_H)-2:
      obs.set_label("V2 (bkg > {:.0f}%)".format(threshold*100))

# Color background of image depending on UHS or SHS
for seg in range(len(UHS_seg)):
   ax1.axvspan(UHS_seg[seg][0], UHS_seg[seg][1], alpha=0.25, color='red')

# Plot simulation results
for file in range(num_data_files_H):
   ax1.scatter(rad2year(data_H[file][:,0]), data_H[file][:,1], s=80, marker=markers[file], c=colors[file], label=labels[file])

# Bottom axis
ax1.set_xlabel('Year', fontsize=20)
ax1.set_ylabel("H$^+$ Rate (s$^{-1}$)", fontsize=20)
for i in range(np.size(V2_path)):
   if V2_year[i] < 2007.0:
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
ax1.set_ylim(0.1,0.5)
ax1.legend(fontsize=20)

# Top axis
ax2.set_xlim(V2_path[idx_left], V2_path[idx_right])
ax2.set_xticks([84.0, 92.0, 100.0, 108.0, 116.0])
ax2.tick_params(labelsize=20)
ax2.set_xlabel("Radial Distance (au)", fontsize=20)

# Vertical lines
y_bot, y_top = ax1.get_ylim()
y_mid = a_bot * y_bot + a_top * y_top
ax1.annotate("TS", (2007.8,y_mid), fontsize=24)
ax1.axvline(2007.67, color='k', linestyle='--', linewidth=2)
ax1.annotate("MAX", (2012.05,y_mid), color='m', fontsize=24)
ax1.axvline(2013.05, color='m', linestyle=':', linewidth=2)

ax1 = fig.add_subplot(313, projection='rectilinear')
ax2 = ax1.twiny()

# Iterate over observation segments based on background threshold
seg_color = "k"
for seg in range(len(cut_idx_e)-1, 0, -1):
   obs, = ax1.semilogy(V2_rate_e[cut_idx_e[seg-1]:cut_idx_e[seg]+1,0], V2_rate_e[cut_idx_e[seg-1]:cut_idx_e[seg]+1,1],
                       color=seg_color, linewidth=2, zorder=0)
   if seg_color == "k":
      seg_color = "c"
   else:
      seg_color = "k"
   if seg == len(cut_idx_e)-1:
      obs.set_label("V2 (bkg < {:.0f}%)".format(threshold*100))
   elif seg == len(cut_idx_e)-2:
      obs.set_label("V2 (bkg > {:.0f}%)".format(threshold*100))

# Color background of image depending on UHS or SHS
for seg in range(len(UHS_seg)):
   ax1.axvspan(UHS_seg[seg][0], UHS_seg[seg][1], alpha=0.25, color='red')

# Plot simulation results
for file in range(num_data_files_e):
   ax1.scatter(rad2year(data_e[file][:,0]), data_e[file][:,1], s=80, marker=markers[file], c=colors[file], label=labels[file])

# Bottom axis
ax1.set_xlabel('Year', fontsize=20)
ax1.set_ylabel("e$^-$ Rate (s$^{-1}$)", fontsize=20)
for i in range(np.size(V2_path)):
   if V2_year[i] < 2007.0:
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
ax1.set_ylim(2.0e-4,2.0e-1)
ax1.legend(fontsize=20)

# Top axis
ax2.set_xlim(V2_path[idx_left], V2_path[idx_right])
ax2.set_xticks([84.0, 92.0, 100.0, 108.0, 116.0])
ax2.tick_params(labelsize=20)
ax2.set_xlabel("Radial Distance (au)", fontsize=20)

# Vertical lines
y_bot, y_top = ax1.get_ylim()
y_mid = np.exp(a_bot * np.log(y_bot) + a_top * np.log(y_top))
ax1.annotate("TS", (2007.8,y_mid), fontsize=24)
ax1.axvline(2007.67, color='k', linestyle='--', linewidth=2)
ax1.annotate("MAX", (2012.05,y_mid), color='m', fontsize=24)
ax1.axvline(2013.05, color='m', linestyle=':', linewidth=2)

# plt.show()
if panel == "left":
   plt.savefig("HS_rate_sim_1-4.png", dpi=300)
else:
   plt.savefig("HS_rate_sim_5-8.png", dpi=300)
plt.close(fig)