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

dat_r25 = np.genfromtxt("r25/norm_realtime_Hop+Onsite+SOC_bias_0.0_dt_0.01_rank_25.txt")
dat_r35 = np.genfromtxt("r35/norm_realtime_Hop+Onsite+SOC_bias_0.0_dt_0.01_rank_35.txt")
dat_r50 = np.genfromtxt("r50/norm_realtime_Hop+Onsite+SOC_bias_0.0_dt_0.01_rank_50.txt")

plt.plot(dat_r25[1::,0],np.abs(1-dat_r25[1::,1]),color="k",linewidth=5,label=r"$r=25$")
plt.plot(dat_r35[1::,0],np.abs(1-dat_r35[1::,1]),color="red",linewidth=5,label=r"$r=35$")
plt.plot(dat_r50[1::,0],np.abs(1-dat_r50[1::,1]),color="blue",linewidth=5,label=r"$r=50$")

#plt.xlim()
plt.yscale("log")
plt.ylim(5*1e-6,3*1e-4)
#plt.yticks([i*1e-6 for i in range(5,51)],labels=[str(i) for i in range(2,9)])
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$\vert 1-\langle\!\langle I\vert \rho\rangle\!\rangle\vert$")

plt.legend(loc="upper left")
plt.savefig("Norm_leads_phonons_spinUncomb_hop+onsite+soc.png")
