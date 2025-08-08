# Code to plot and fit LISM spectra

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import quad
import sys

# Broken power law
def log_broken_pow_law(T, J, T_b, a1, a2, d):
   return np.log(J * (T / 1000)**a1 / (1.0 + (T / T_b)**(a2/d))**d)

def broken_pow_law(T, J, T_b, a1, a2, d):
   return np.exp(log_broken_pow_law(T, J, T_b, a1, a2, d))

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
if specie_label == "He":
   voyager2_2018 = np.loadtxt("data/V2_He_flux_2018_001_2018_305.dat")
   voyager2_2019 = np.loadtxt("data/V2_He_flux_2019_070_2019_159.dat")
   bess = np.loadtxt("data/BESS_He_flux_1997_2002.dat")
   energy_V = voyager2_2018[:,0]
   flux_V = voyager2_2018[:,1]
   label_V = "V2 (2018)"
   energy_W = voyager2_2019[:,0]
   flux_W = voyager2_2019[:,1]
   label_W = "V2 (2019)"
   energy_X = bess[:,0] * 1000
   flux_X = bess[:,1] / 1000
   label_X = "BESS (1997-2002)"
   low_energy = energy_W[0]
   high_energy = energy_X[-1]
   guess_params = [100.0, 60.0, 0.5, 3.0, 5.0]
   low_params = [1.0, 1.0, 0.0, 2.0, 1.0]
   high_params = [200.0, 1000.0, 1.0, 4.0, 10.0]
   energy_label = "MeV/nuc"
elif specie_label == "e":
   voyager1 = np.loadtxt("data/V1_e_flux_2012_342_2015_181.dat")
   ams = np.loadtxt("data/AMS_e_flux_2011_139_2013_330.dat")
   energy_V = voyager1[:,0]
   flux_V = voyager1[:,1]
   label_V = "V1 (2012-2015)"
   energy_X = ams[:,0] * 1000
   flux_X = ams[:,1] / 1000
   label_X = "AMS (2011-2013)"
   low_energy = energy_V[0]
   high_energy = energy_X[-1]
   guess_params = [10.0, 100.0, -2.0, 3.0, 2.0]
   low_params = [0.1, 1.0, -3.0, 1.0, 1.0]
   high_params = [100.0, 1000.0, 0.0, 4.0, 5.0]
   energy_label = "MeV"

# Fit data
energy_F = np.exp(np.linspace(np.log(low_energy), np.log(high_energy), num=100))
energy = np.concatenate((energy_V, energy_X))
flux = np.concatenate((flux_V, flux_X))
opt_params, cov = curve_fit(log_broken_pow_law, energy, np.log(flux),
                            p0=guess_params, maxfev=10000,
                            bounds=(low_params, high_params))

# Import data from results
perc = 50
mod_spec = np.loadtxt("../results/HS_mod_spec_{:s}/HS_mod_parker_{:d}_pct_spec_comp.dat".format(specie_label,perc))
mod_spec[:,0] = mod_spec[:,0] * 1000 # GeV -> MeV

# Plot
fig = plt.figure(figsize=(8, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')
energy = np.exp(np.linspace(np.log(1.0), np.log(1000.0), num=100))
ax.loglog(energy, broken_pow_law(energy, *opt_params), color="tab:red", linewidth=2)
ax.loglog(mod_spec[:,0], mod_spec[:,1], color="tab:blue", label="Unmodulated", linewidth=3, linestyle='--')
ax.loglog(mod_spec[:,0], mod_spec[:,2], color="tab:orange", label="Modulated", linewidth=3)
ax.set_xlabel("Kinetic Energy ({:s})".format(energy_label), fontsize = 20)
ax.set_ylabel(specie_label + " Flux (m$^2$ s sr " + energy_label + ")$^{-1}$", fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize = 20)

if specie_label == "He":
   E1 = 130.0
   E2 = 460.0
elif specie_label == "e":
   E1 = 3.0
   E2 = 100.0
ax.axvline(E1, color='k', linestyle='--', linewidth=2)
ax.axvline(E2, color='k', linestyle='--', linewidth=2)

plt.show()
