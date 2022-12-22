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

dat1 = np.genfromtxt("L1_gPhU0.0_gPhD0.0_V0.0_w00.1_X0_1.0_dt0.002_1_phonons.dat")
t = dat1[:,0]

w = 0.1

X0 = 1.0
P0 = 0.0
A = X0/2 + P0/(2*1j*w)
B = -1/(2*1j*w)*(P0-1j*w*X0)

X_analyt = A*np.exp(1j*w*t)+B*np.exp(-1j*w*t)
P_analyt = 1j*w*(A*np.exp(1j*w*t)-B*np.exp(-1j*w*t))

ax1.plot(t,dat1[:,1],color="black",linewidth=6,label=r"$X$ numerical")
ax1.plot(t,X_analyt,color="red",linewidth=4,linestyle="dashed",label=r"$X$ analytical")
ax1.plot(t,dat1[:,2],color="blue",linewidth=6,label=r"$P$ numerical")
ax1.plot(t,P_analyt,color="orange",linewidth=4,linestyle="dashed",label=r"$P$ analytical")

ax1.legend(loc="lower right")
ax1.set_ylim(-1.1,1.1)

ax1.set_xlabel(r"$t$",labelpad=-3)
plt.savefig("Decoupled_vibrational_dynamics.pdf")
