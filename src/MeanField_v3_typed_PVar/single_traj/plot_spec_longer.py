import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

plt.style.use("classic")
plt.rc("text", usetex=True)
plt.rc("font",family="serif",serif="Times")
fs = 21
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*0.2*4.135668*1e3)
fac_e = 0.002*0.2

plt.rc("axes",labelsize=fs)
plt.rc("xtick",labelsize=fs-2)
plt.rc("ytick",labelsize=fs-2)

px = 0.14
py = 0.155
sx = 0.85
sy = 0.78
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(8,4.5))
ax = plt.axes([px,py,sx,sy])
sp = "P"
dat_sp = np.genfromtxt("spec_P_L3x7.dat")
dat_el = np.genfromtxt("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons_es.dat")
N = len(dat_sp[0,1::])
x = dat_sp[0::10,0]*fac_t
mn = dat_el[:,1::].min()
mx = dat_el[:,1::].max()
norm = plt.Normalize(mn, mx)
for i in range(N):
    y = dat_sp[0::10,i+1]*fac_e
    c = dat_el[0::10,i+1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='viridis', norm=norm)
    # Set the values used for colormapping
    lc.set_array(c)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)


cbar = fig.colorbar(line, ax=ax)
cbar.set_label("adiabatic state population",labelpad=10)
cbar.ax.tick_params(labelsize=fs-2)
ax.set_xlim(segments[:,:,0].min(), segments[:,:,0].max())
ax.set_ylim(-0.5,0.0)
ax.set_ylabel("energy "+r"$E\,[\mathrm{eV}]$")
ax.set_xlabel("time "+r"$t\, [\mathrm{ps}]$")
ax.annotate(r"$A=+\frac{1}{2}$",(0.05,0.9),xycoords="axes fraction",fontsize=fs)
plt.savefig(f"L3x6_elec_es_spec_{sp}_singleTraj_longer.pdf")
