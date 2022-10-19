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

dat_r25 = np.genfromtxt("r60/n_site_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_60.txt")[1::,:]
dat_r35 = np.genfromtxt("r90/n_site_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_90.txt")[1::,:]
dat_r50 = np.genfromtxt("r120/n_site_realtime_Hop+Onsite+Holstein_bias_0.0_dt_0.01_rank_120.txt")[1::,:]
dat_iT  = np.genfromtxt("../../tdvp/N3_t01.00_tso0.000_tPh0.100_dt0.00500.dat")[1::2,1]



d1 = np.abs(dat_iT-dat_r25[:,1])
d2 = np.abs(dat_iT-dat_r35[:,1])
d3 = np.abs(dat_iT-dat_r50[:,1])


plt.plot(dat_r25[:,0],d1,color="k",linewidth=5,label=r"$r=25$")
plt.plot(dat_r35[:,0],d2,color="red",linewidth=5,label=r"$r=80$")
plt.plot(dat_r50[:,0],d3,color="blue",linewidth=5,linestyle="dashed",label=r"$r=120$")


#plt.xlim()
plt.xscale("log")
plt.yscale("log")
#plt.ylim(2*1e-6,0.81e-2)
#plt.yticks([i*1e-5 for i in range(2,9)],labels=[str(i) for i in range(2,9)])
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$\vert n_{1\uparrow}^{\mathrm{GSE-TDVP}}-n_{1\uparrow}^{\mathrm{1TDVP}}\vert$")

#plt.legend(loc="upper left")
plt.savefig("Compare_leads_phonons_spinUncomb_hop+onsite+holstein.pdf")
