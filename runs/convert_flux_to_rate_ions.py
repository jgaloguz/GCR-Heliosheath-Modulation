# Code to convert flux to rate for ion species in V2 measurements

# Import library
import numpy as np

# Function to import data
def ImportData(filename, n_channels):
   labels = []
   year_ions = []
   fluxes = [[] for _ in range(n_channels)]
   errors = [[] for _ in range(n_channels)]
   counts = [[] for _ in range(n_channels)]
   file = open(filename, "r")
# Header line
   file.readline()
# Labels
   for l in range(n_channels):
      line = file.readline()
      labels.append(line)
# Data
   line = file.readline()
   while line:
      data_str = line.split()
      year_ions.append(float(data_str[0]))
      for c in range(n_channels):
         fluxes[c].append(float(data_str[1+3*c]))
         errors[c].append(float(data_str[1+3*c+1]))
         counts[c].append(float(data_str[1+3*c+2]))
      line = file.readline()
   file.close()
   return labels, fluxes, errors, counts, year_ions

def Integrate(labels, fluxes, errors, n_channels, geom_fact):
   flux = np.zeros(np.size(fluxes,1))
   ferr = np.zeros(np.size(errors,1))
   for c in range(n_channels):
      flux = flux + np.array(fluxes[c]) * (float(labels[c][11:18]) - float(labels[c][1:8]))
      ferr = ferr + np.array(errors[c]) * (float(labels[c][11:18]) - float(labels[c][1:8]))
   return flux * geom_fact, ferr * geom_fact

# Import and format data
n_chan_H = 7
labels_H, fluxes_H, errors_H, counts_H, year_H = ImportData("data/V2_H_flux.dat", n_chan_H)
n_chan_He = 6
labels_He, fluxes_He, errors_He, counts_He, year_He = ImportData("data/V2_He_flux.dat", n_chan_He)

# Add up contribution from all channels
rate_H, rerr_H = Integrate(labels_H, fluxes_H, errors_H, n_chan_H, 1.678e-4)
rate_He, rerr_He = Integrate(labels_He, fluxes_He, errors_He, n_chan_He, 1.362e-4)

# Save arrays to file
np.savetxt("data/V2_H_rate.dat", np.vstack((year_H, rate_H, rerr_H)).T)
np.savetxt("data/V2_He_rate.dat", np.vstack((year_He, rate_He, rerr_He)).T)
