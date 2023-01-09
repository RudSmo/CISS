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
ax1 = plt.axes([0.1,0.147,0.86,0.82])
#ax1 = plt.axes([0.1,0.5,0.86,0.45])
#ax2 = plt.axes([0.1,0.12,0.86,0.45])

dat1 = np.genfromtxt("fS3_nL1_tso0.001_P0.5_currents.dat")
dat2 = np.genfromtxt("fS3_nL1_tso0.001_P0.0_currents.dat")
dat3 = np.genfromtxt("fS3_nL1_tso0.001_P-0.5_currents.dat")


t0 = 0.04
Gamma = 0.25*t0

ax1.plot(dat1[:,0],dat1[:,1]/Gamma,color="orange",linewidth=6,label=r"$P=0.5$")
ax1.plot(dat2[:,0],dat2[:,1]/Gamma,color="black",linewidth=5,linestyle="dashed",dashes=(9,10),label=r"$P=0.0$")
ax1.plot(dat3[:,0],dat3[:,1]/Gamma,color="blue",linewidth=4,linestyle="dashed",dashes=(10,10),label=r"$P=-0.5$")

ax1.legend(loc="lower right")

ax1.set_xlim(-0.1*t0,35*t0)
ax1.set_ylabel(r"$I \, [\Gamma_0]$")
#ax1.tick_params(bottom=False)
ax1.set_xlabel(r"$V \, [\mathrm{V}]$")

#ax2.plot(dat1[:,0],dat1[:,1]-dat3[:,1],color="k",linewidth=5)
#ax2.set_xlabel(r"$V \, [t_0]$")
#ax2.set_ylabel(r"$I(+)-I(-)$")
plt.savefig("LB_1x3.pdf")
