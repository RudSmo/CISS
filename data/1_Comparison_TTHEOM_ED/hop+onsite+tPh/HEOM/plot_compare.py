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

dat1 = np.genfromtxt("Leads_Phonons_spinUncomb/r120/n_site_realtime_Hop+Onsite+tPh_bias_0.0_dt_0.01_rank_120.txt")
dat2 = np.genfromtxt("Leads_Phonons_spinUncomb_dummyInds/r120/n_site_realtime_Hop+Onsite_bias_0.0_dt_0.01_rank_120.txt")
dat3 = np.genfromtxt("Leads_PhononsAlt_spinUncomb/r128/n_site_realtime_Hop+Onsite_bias_0.0_dt_0.01_rank_128.txt")
dat4 = np.genfromtxt("Leads_PhononsAlt_spinComb/r128/n_site_up_realtime_bias_0.0_dt_0.01_rank_128.txt")

plt.plot(dat1[1::,0],dat1[1::,1],color="k",linewidth=7,label="		")
plt.plot(dat2[1::,0],dat2[1::,1],color="red",linestyle="dashed",linewidth=6,label="		")
plt.plot(dat3[1::,0],dat3[1::,1],color="blue",linestyle="dashed",dashes=(7,10),linewidth=5,label=r"$\qquad \qquad$")
plt.plot(dat4[1::,0],dat4[1::,1],color="green",linewidth=5,label=r"$\qquad \qquad$")

#plt.xlim()
#plt.yscale("log")
plt.ylim(-0.1,1.2)
#plt.yticks([i*1e-6 for i in range(5,51)],labels=[str(i) for i in range(2,9)])
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$n_{1\uparrow}$")

#plt.legend(loc="upper right",ncol=2,columnspacing=5,borderpad=0.4)
plt.savefig("Compare_MPS_toplogy_hop+onsite+tPh.png")
