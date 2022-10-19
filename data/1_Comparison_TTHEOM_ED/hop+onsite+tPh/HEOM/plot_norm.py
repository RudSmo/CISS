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


fig = plt.figure(figsize=(10,8))

dat1 = np.genfromtxt("Leads_Phonons_spinUncomb/r120/norm_realtime_Hop+Onsite+tPh_bias_0.0_dt_0.01_rank_120.txt")
dat2 = np.genfromtxt("Leads_Phonons_spinUncomb_dummyInds/r120/norm_realtime_Hop+Onsite_bias_0.0_dt_0.01_rank_120.txt")

plt.plot(dat1[1::,0],np.abs(1-dat1[1::,1]),color="k",linewidth=7,label="		")
plt.plot(dat2[1::,0],np.abs(1-dat2[1::,1]),color="red",linestyle="dashed",linewidth=6,label="		")

#plt.xlim()
plt.yscale("log")
#plt.ylim(-0.1,1.2)
#plt.yticks([i*1e-6 for i in range(5,51)],labels=[str(i) for i in range(2,9)])
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$\vert 1-\langle\!\langle I\vert \rho\rangle\!\rangle\vert$")

#plt.legend(loc="lower right",ncol=2,columnspacing=5,borderpad=0.8)
plt.savefig("Norm_MPS_toplogy_hop+onsite+tPh.png")
