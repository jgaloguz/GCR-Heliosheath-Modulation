# Code to plot diffusion calculations

# Import modules
import numpy as np
import matplotlib.pyplot as plt
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

# Import data
kappa_SIM_rig = np.loadtxt("../results/kappa_SIM_rig_{:s}.dat".format(specie_label))
kappa_SIM_V2 = np.loadtxt("../results/kappa_SIM_V2_{:s}.dat".format(specie_label))
kappa_SOQLT_rig_SHS = np.loadtxt("../results/kappa_SOQLT_rig_{:s}_SHS.dat".format(specie_label))
kappa_SOQLT_V2_SHS = np.loadtxt("../results/kappa_SOQLT_V2_{:s}_SHS.dat".format(specie_label))
kappa_UNLT_rig_SHS = np.loadtxt("../results/kappa_UNLT_rig_{:s}_SHS.dat".format(specie_label))
kappa_UNLT_V2_SHS = np.loadtxt("../results/kappa_UNLT_V2_{:s}_SHS.dat".format(specie_label))
kappa_SOQLT_rig_UHS = np.loadtxt("../results/kappa_SOQLT_rig_{:s}_UHS.dat".format(specie_label))
kappa_SOQLT_V2_UHS = np.loadtxt("../results/kappa_SOQLT_V2_{:s}_UHS.dat".format(specie_label))
kappa_UNLT_rig_UHS = np.loadtxt("../results/kappa_UNLT_rig_{:s}_UHS.dat".format(specie_label))
kappa_UNLT_V2_UHS = np.loadtxt("../results/kappa_UNLT_V2_{:s}_UHS.dat".format(specie_label))

# Plot
fig = plt.figure(figsize=(15, 10), layout='tight')
ax1 = fig.add_subplot(211, projection='rectilinear')

ax1.loglog(kappa_SOQLT_rig_SHS[:,0], kappa_SOQLT_rig_SHS[:,1], color = 'tab:red',
           linestyle="", marker="s", markersize=6, label='SOQLT')
ax1.loglog(kappa_SIM_rig[:,0], kappa_SIM_rig[:,1], color = 'tab:blue',
           linestyle="-", linewidth = 2, label='EMP $\\parallel$')
ax1.loglog(kappa_UNLT_rig_SHS[:,0], kappa_UNLT_rig_SHS[:,1], color = 'tab:orange',
           linestyle="", marker="s", markersize=6, label='UNLT')
ax1.loglog(kappa_SIM_rig[:,0], kappa_SIM_rig[:,3], color = 'tab:cyan',
           linestyle="-", linewidth = 2, label='EMP $\\perp$')
ax1.set_xlabel('R (GV)', fontsize=20)
ax1.set_ylabel('$\\kappa_{\\parallel, \\perp}$ (cm$^2$ s$^{-1}$)', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
ax1.legend(fontsize=20)

ax2 = fig.add_subplot(212, projection='rectilinear')

ax2.semilogy(kappa_SOQLT_V2_SHS[:,0], kappa_SOQLT_V2_SHS[:,1], color = 'tab:red',
             linestyle="", marker="s", markersize=3, label='SOQLT')
ax2.semilogy(kappa_SOQLT_V2_UHS[:,0], kappa_SOQLT_V2_UHS[:,1], color = 'tab:red',
             linestyle="", marker="s", markersize=3)
ax2.semilogy(kappa_SIM_V2[:,0], kappa_SIM_V2[:,1], color = 'tab:blue',
             linestyle="-", linewidth = 2, label='EMP $\\parallel$')
ax2.semilogy(kappa_UNLT_V2_SHS[:,0], kappa_UNLT_V2_SHS[:,1], color = 'tab:orange',
             linestyle="", marker="s", markersize=3, label='UNLT')
ax2.semilogy(kappa_UNLT_V2_UHS[:,0], kappa_UNLT_V2_UHS[:,1], color = 'tab:orange',
             linestyle="", marker="s", markersize=3)
ax2.semilogy(kappa_SIM_V2[:,0], kappa_SIM_V2[:,2], color = 'tab:cyan',
             linestyle="-", linewidth = 2, label='EMP $\\perp$')
ax2.set_xlabel('r (au)', fontsize=20)
ax2.set_ylabel('$\\kappa_{\\parallel, \\perp}$ (cm$^2$ s$^{-1}$)', fontsize=20)
ax2.set_xlim(80.0,120.0)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.axvline(83.6, color='k', linestyle='--', linewidth=2)
ax2.annotate("TS", (84.0, 3.0e22), fontsize=24)
ax2.axvline(119.0, color='k', linestyle='--', linewidth=2)
ax2.annotate("HP", (117.0, 3.0e22), fontsize=24)

plt.show()
plt.close(fig)