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
   return np.log(J * (T / 1000)**a1 / (1.0 + (T / T_b)**((a1-a2)/d))**d)

def broken_pow_law(T, J, T_b, a1, a2, d):
   return np.exp(log_broken_pow_law(T, J, T_b, a1, a2, d))

# Import data
He_voyager2 = np.loadtxt("data/V2_He_flux_2019_070_2019_159.dat")
He_bess = np.loadtxt("data/BESS_He_flux_1997_2002.dat")
He_energy_V = He_voyager2[:,0]
He_flux_V = He_voyager2[:,1]
He_energy_X = He_bess[:,0] * 1000
He_flux_X = He_bess[:,1] / 1000
He_low_energy = He_energy_V[0]
He_high_energy = He_energy_X[-1]
He_guess_params = [5.0,   60.0,   0.5, -3.0, 5.0]
He_low_params =   [1.0,   1.0,    0.0, -4.0, 1.0]
He_high_params =  [100.0, 1000.0, 1.0, -2.0, 10.0]
H_voyager2 = np.loadtxt("data/V2_H_flux_2019_070_2019_159.dat")
H_bess = np.loadtxt("data/BESS_H_flux_1997_2002.dat")
H_energy_V = H_voyager2[:,0]
H_flux_V = H_voyager2[:,1]
H_energy_X = H_bess[:,0] * 1000
H_flux_X = H_bess[:,1] / 1000
H_low_energy = H_energy_V[0]
H_high_energy = H_energy_X[-1]
H_guess_params = [50.0,   60.0,   0.5, -3.0, 5.0]
H_low_params =   [10.0,   1.0,    0.0, -4.0, 1.0]
H_high_params =  [1000.0, 1000.0, 2.0, -2.0, 10.0]
H_energy_label = "MeV"
e_voyager1 = np.loadtxt("data/V1_e_flux_2012_342_2015_181.dat")
e_ams = np.loadtxt("data/AMS_e_flux_2011_139_2013_330.dat")
e_energy_V = e_voyager1[:,0]
e_flux_V = e_voyager1[:,1]
e_energy_X = e_ams[:,0] * 1000
e_flux_X = e_ams[:,1] / 1000
e_low_energy = e_energy_V[0]
e_high_energy = e_energy_X[-1]
e_guess_params = [10.0,  100.0,  -2.0, -3.0, 2.0]
e_low_params =   [0.1,   1.0,    -3.0, -4.0, 1.0]
e_high_params =  [100.0, 1000.0, 0.0,  -2.0, 5.0]

# Fit data
He_energy_F = np.exp(np.linspace(np.log(He_low_energy), np.log(He_high_energy), num=100))
He_energy = np.concatenate((He_energy_V, He_energy_X))
He_flux = np.concatenate((He_flux_V, He_flux_X))
He_opt_params, cov = curve_fit(log_broken_pow_law, He_energy, np.log(He_flux),
                               p0=He_guess_params, maxfev=10000,
                               bounds=(He_low_params, He_high_params))

H_energy_F = np.exp(np.linspace(np.log(H_low_energy), np.log(H_high_energy), num=100))
H_energy = np.concatenate((H_energy_V, H_energy_X))
H_flux = np.concatenate((H_flux_V, H_flux_X))
H_opt_params, cov = curve_fit(log_broken_pow_law, H_energy, np.log(H_flux),
                              p0=H_guess_params, maxfev=10000,
                              bounds=(H_low_params, H_high_params))

e_energy_F = np.exp(np.linspace(np.log(e_low_energy), np.log(e_high_energy), num=100))
e_energy = np.concatenate((e_energy_V, e_energy_X))
e_flux = np.concatenate((e_flux_V, e_flux_X))
e_opt_params, cov = curve_fit(log_broken_pow_law, e_energy, np.log(e_flux),
                              p0=e_guess_params, maxfev=10000,
                              bounds=(e_low_params, e_high_params))

# Plot
fig = plt.figure(figsize=(8, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.loglog(He_energy_V, He_flux_V, linestyle="", marker="^",
          color="tab:blue", markersize=8)
ax.loglog(He_energy_X, He_flux_X, linestyle="", marker="^",
          color="tab:green", markersize=8)
ax.loglog(He_energy_F, broken_pow_law(He_energy_F, *He_opt_params),
          color="tab:red", label="He$^{2\\!+}$", linewidth=2, linestyle="-")
ax.loglog(H_energy_V, H_flux_V, linestyle="", marker="s",
          color="tab:blue", markersize=8)
ax.loglog(H_energy_X, H_flux_X, linestyle="", marker="s",
          color="tab:green", markersize=8)
ax.loglog(H_energy_F, broken_pow_law(H_energy_F, *H_opt_params),
          color="tab:red", label="H$^{+}$", linewidth=2, linestyle="--")
ax.loglog(e_energy_V, e_flux_V, linestyle="", marker="o",
          color="tab:blue", markersize=8)
ax.loglog(e_energy_X, e_flux_X, linestyle="", marker="o",
          color="tab:green", markersize=8)
ax.loglog(e_energy_F, broken_pow_law(e_energy_F, *e_opt_params),
          color="tab:red", label="e$^{-}$", linewidth=2, linestyle=":")

ax.set_xlabel("Kinetic Energy (MeV or MeV/nuc)", fontsize = 20)
ax.set_ylabel("Flux (m$^2$ s sr MeV or MeV/nuc)$^{-1}$", fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize = 20)

plt.show()
plt.close(fig)

# Compute and report integrated fluxes
print("\nOptimal Fit Parameters He:")
print(He_opt_params)
E1 = 130.0
E2 = 460.0
geom_factor = 1.50e-4 # m^2 s (Cummings et al 2016)
I = quad(broken_pow_law, E1, E2, args=tuple(He_opt_params))
S = I[0] * geom_factor
print("Integrated flux over {:.0f}-{:.0f} {:s} = {:.3f} s^-1".format(E1, E2, "MeV/nuc", S))

print("\nOptimal Fit Parameters H:")
print(H_opt_params)
E1 = 130.0
E2 = 345.0
geom_factor = 1.68e-4 # m^2 s  (Cummings et al 2016)
I = quad(broken_pow_law, E1, E2, args=tuple(H_opt_params))
S = I[0] * geom_factor
print("Integrated flux over {:.0f}-{:.0f} {:s} = {:.3f} s^-1".format(E1, E2, "MeV", S))

print("\nOptimal Fit Parameters e:")
print(e_opt_params)
N = 100
E1 = 3.0
E2 = 100.0
energy_array = np.exp(np.linspace(np.log(E1), np.log(E2), num=N))
response_datafile = np.loadtxt("data/TET_D123_response_m2sr.dat")
response_func = interp1d(response_datafile[:,0], response_datafile[:,1])
S = 0.0
for i in range(N-1):
   E = 0.5 * (energy_array[i] + energy_array[i+1])
   dE = (energy_array[i+1] - energy_array[i])
   S += response_func(E) * broken_pow_law(E, *e_opt_params) * dE
print("Integrated flux over {:.0f}-{:.0f} {:s} = {:.3f} s^-1".format(E1, E2, "MeV", S))

