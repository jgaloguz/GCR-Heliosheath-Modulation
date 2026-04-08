# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Arrow parameters
X0 = -0.6
Y0 = -1.5
N = 4
dang = 0.5 * np.pi / (N - 1)

# x/o parameters
M = 4
X_x = np.linspace(-1.2, -0.2, num=4)
Y_x = np.linspace(-0.2,  0.8, num=4)
XX_x, YY_x = np.meshgrid(X_x, Y_x)

X_o = np.linspace( 0.2,  1.2, num=4)
Y_o = np.linspace(-0.2,  0.8, num=4)
XX_o, YY_o = np.meshgrid(X_o, Y_o)

# Draw
fig = plt.figure(figsize=(8, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

xs = ax.scatter(XX_x, YY_x, s=100.0, c="tab:green", marker="x", label="$B_z$ (< 0)")
xo = ax.scatter(XX_o, YY_o, s=100.0, c="tab:green", marker="o", label="$B_z$ (> 0)")

n_cross = 0
label_cross = True
label_no_cross = True
for ang in range(N):
   R = 1.0
   T = ang * dang
   DX = R * np.cos(T)
   DY = R * np.sin(T)
   if X0+DX > 0.0:
      n_cross = n_cross + 1
      c = "tab:blue"
      if label_cross:
         lab="gradient\nacross sheet"
         label_cross = False
      else:
         lab=""
   else:
      c = "tab:red"
      if label_no_cross:
         lab="gradient\nbeside sheet"
         label_no_cross = False
      else:
         lab=""
   ag = ax.arrow(X0, Y0, DX, DY, color=c,
                 width=0.02, length_includes_head=True, label=lab)
ad = ax.arrow(X0, Y0, 0.0, -n_cross/N, color="tab:purple",
              width=0.02, length_includes_head=True, label="drift")
pp = ax.scatter(X0, Y0, s=100.0, c="k", marker="o", label="particle")

ax.axvline(0.0, color='k', linestyle="--", linewidth=5, label="sheet")
ax.set_xlabel("x", fontsize = 20)
ax.set_ylabel("y", fontsize = 20)
ax.set_xlim(-1.5,1.5)
ax.set_ylim(-2.3,1.2)
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])

ax.annotate("$B_z$ (< 0)", (-0.9, 1.0), color="tab:green", fontsize=20)
ax.annotate("$B_z$ (> 0)", ( 0.5, 1.0), color="tab:green", fontsize=20)
ax.annotate("grad\nacross\nsheet", (X0+R+0.1, Y0-0.17), color="tab:blue", fontsize=20)
ax.annotate("grad\nalong\nsheet", (X0-0.45, Y0+0.3), color="tab:red", fontsize=20)
ax.annotate("drift", (X0-0.3, Y0-n_cross/N+0.02), color="tab:purple", fontsize=20)
ax.annotate("particle", (X0-0.5, Y0-0.03), color="k", fontsize=20)
ax.annotate("current sheet", (0.1, -2.1), color="k", fontsize=20)


# ax.legend(loc=4, fontsize=20)

plt.savefig("test.png")
plt.show()
plt.close(fig)