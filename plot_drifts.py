# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# which quantity to plot
which_plot = 0

# Import plotting parameters from file
drift_params_file = open("params_drifts.txt", "r")
Nx = int(drift_params_file.readline())
Nz = int(drift_params_file.readline())
X0 = float(drift_params_file.readline())
Z0 = float(drift_params_file.readline())
DX = float(drift_params_file.readline())
DZ = float(drift_params_file.readline())
drift_params_file.close()

X = np.linspace(X0, X0+DX, num=Nx+1)
Z = np.linspace(Z0, Z0+DZ, num=Nz+1)
XX, ZZ = np.meshgrid(X, Z)
dx = 0.5*DX/Nx
dz = 0.5*DZ/Nz
X_q = np.linspace(X0+dx, X0+DX-dx, num=Nx)
Z_q = np.linspace(Z0+dz, Z0+DZ-dz, num=Nz)
XX_q, ZZ_q = np.meshgrid(X_q, Z_q)
VX = np.zeros((Nz,Nx))
VY = np.zeros((Nz,Nx))
VZ = np.zeros((Nz,Nx))
VM = np.zeros((Nz,Nx))
RL = np.zeros((Nz,Nx))
KX = np.zeros((Nz,Nx))
KY = np.zeros((Nz,Nx))
KZ = np.zeros((Nz,Nx))
KM = np.zeros((Nz,Nx))

geometry = "flat"
file = open("../results/gcr_drifts.dat", 'r')
for j in range(Nx):
   for i in range(Nz):
      line = file.readline().split()
      VX[i,j] = float(line[0])
      VY[i,j] = float(line[1])
      VZ[i,j] = float(line[2])
      VM[i,j] = np.sqrt(VX[i,j]**2 + VY[i,j]**2 + VZ[i,j]**2)
      RL[i,j] = float(line[3])
      if which_plot != 0:
         KX[i,j] = float(line[4])
         KY[i,j] = float(line[5])
         KZ[i,j] = float(line[6])
         KM[i,j] = np.sqrt(KX[i,j]**2 + KY[i,j]**2 + KZ[i,j]**2)
      
file.close()

# Plots
if which_plot == 0:
   title = "v_d"
   UM = VM
   UX = VX
   UY = VY
   UZ = VZ
   pic_name = "drift"
else:
   title = "\\nabla \\cdot \\kappa"
   UM = KM
   UX = KX
   UY = KY
   UZ = KZ
   pic_name = "divK"

fig = plt.figure(figsize=(14, 12), layout='tight')

ax1 = fig.add_subplot(221, projection='rectilinear')
ax1.set_xlabel("AU", fontsize = 16)
ax1.set_ylabel("AU", fontsize = 16)
ax1.set_title("$|"+ title +"|$", fontsize = 16)
ax1.tick_params(labelsize=16)
hm = ax1.pcolormesh(XX, ZZ, UM)
cb1 = fig.colorbar(hm, ax=ax1)
cb1.ax.tick_params(labelsize=16)

ax2 = fig.add_subplot(222, projection='rectilinear')
ax2.set_xlabel("AU", fontsize = 16)
ax2.set_ylabel("AU", fontsize = 16)
ax2.set_title("$(" + title + ")_x$", fontsize = 16)
ax2.tick_params(labelsize=16)
hm = ax2.pcolormesh(XX, ZZ, UX)
cb2 = fig.colorbar(hm, ax=ax2)
cb2.ax.tick_params(labelsize=16)

ax3 = fig.add_subplot(223, projection='rectilinear')
ax3.set_xlabel("AU", fontsize = 16)
ax3.set_ylabel("AU", fontsize = 16)
ax3.set_title("$(" + title + ")_y$", fontsize = 16)
ax3.tick_params(labelsize=16)
hm = ax3.pcolormesh(XX, ZZ, UY)
cb3 = fig.colorbar(hm, ax=ax3)
cb3.ax.tick_params(labelsize=16)

ax4 = fig.add_subplot(224, projection='rectilinear')
ax4.set_xlabel("AU", fontsize = 16)
ax4.set_ylabel("AU", fontsize = 16)
ax4.set_title("$(" + title + ")_z$", fontsize = 16)
ax4.tick_params(labelsize=16)
hm = ax4.pcolormesh(XX, ZZ, UZ)
cb4 = fig.colorbar(hm, ax=ax4)
cb4.ax.tick_params(labelsize=16)

plt.savefig("../results/gcr_drifts.png", dpi=300)
plt.show()
plt.close(fig)

# fig = plt.figure(figsize=(10, 8), layout='tight')
# ax = fig.add_subplot(111, projection='rectilinear')

# qp = ax.streamplot(XX_q, ZZ_q, VX, VZ, density=2)
# ax.set_xlabel("AU")
# ax.set_ylabel("AU")

# plt.show()
# plt.close(fig)