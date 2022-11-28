import numpy as np
import matplotlib.pyplot as plt
import os,fnmatch
from matplotlib import rc

plt.style.use("classic")

ls = 23
lw = 3
rc("font",family="serif",serif="Times")
rc("text",usetex=True)
rc("xtick",labelsize=ls)
rc("ytick",labelsize=ls)
rc("axes",labelsize=ls)
rc("legend",fontsize=ls-5)


fig = plt.figure(figsize=(10,5))
ax1 = plt.axes([0.1,0.152,0.86,0.83])
#ax1 = plt.axes([0.1,0.5,0.86,0.45])
#ax2 = plt.axes([0.1,0.12,0.86,0.45])

dat1 = np.genfromtxt("L1_gPhU0.0_gPhD0.0_dt0.001_1_phonons.dat")
Xrng = dat1[:,0]

omega = 0.1

X_analyt = np.cos(omega**(3/2)*Xrng)
P_analyt = -np.sqrt(omega)*np.sin(omega**(3/2)*Xrng)

ax1.plot(dat1[:,0],dat1[:,1],color="black",linewidth=6,label=r"$X$ numerical")
ax1.plot(Xrng,X_analyt,color="red",linewidth=4,linestyle="dashed",label=r"$X$ analytical")
ax1.plot(dat1[:,0],dat1[:,2],color="blue",linewidth=6,label=r"$P$ numerical")
ax1.plot(Xrng,P_analyt,color="orange",linewidth=4,linestyle="dashed",label=r"$P$ analytical")

ax1.legend(loc="lower right")
ax1.set_ylim(-1.1,1.1)

ax1.set_xlabel(r"$t \, [t_0^{-1}]$",labelpad=-3)
plt.savefig("Decoupled_vibrational_dynamics.pdf")
