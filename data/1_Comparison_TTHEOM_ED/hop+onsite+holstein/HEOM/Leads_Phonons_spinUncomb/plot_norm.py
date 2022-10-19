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
rc("legend",fontsize=ls)


fig = plt.figure(figsize=(10,4))

dat_r25 = np.genfromtxt("r60/norm_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_60.txt")
dat_r35 = np.genfromtxt("r90/norm_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_90.txt")
dat_r50 = np.genfromtxt("r120/norm_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_120.txt")

plt.plot(dat_r25[1::,0],np.abs(1-dat_r25[1::,1]),color="k",linewidth=7,label=r"$r=60$")
plt.plot(dat_r35[1::,0],np.abs(1-dat_r35[1::,1]),color="red",linewidth=4,label=r"$r=90$")
plt.plot(dat_r50[1::,0],np.abs(1-dat_r50[1::,1]),color="blue",linewidth=4,label=r"$r=120$")

#plt.xlim()
plt.yscale("log")
#plt.ylim(2*1e-6,0.81e-2)
#plt.yticks([i*1e-5 for i in range(2,9)],labels=[str(i) for i in range(2,9)])
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$\vert 1-\langle\!\langle I\vert \rho\rangle\!\rangle\vert$",labelpad=5)

plt.legend()
plt.savefig("Norm_leads_phonons_spinUncomb_hop+onsite+holstein.png")
