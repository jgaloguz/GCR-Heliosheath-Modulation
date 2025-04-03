import numpy as np
import matplotlib.pyplot as plt

# Import data
data = np.loadtxt("data/TET_D123_response_m2sr.dat")
energy = data[:,0]
response = data[:,1]

# Plot
fig = plt.figure(figsize=(12, 8), layout='tight')
ax1 = fig.add_subplot(111, projection='rectilinear')

ax1.loglog(energy, response, 'b-', linewidth=3)
ax1.set_xlabel('Energy (MeV)', fontsize=20)
ax1.set_ylabel('D1-D3 Response (m$^2$ sr)', fontsize=20)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)

plt.show()
plt.close(fig)