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


fig = plt.figure(figsize=(10,7))
ax0 = plt.axes([0.13,0.12,0.84,0.85])
ax1 = plt.axes([0.45,0.2,0.5,0.3])

dat_HE1 = np.genfromtxt("TT-HEOM/nd2/np10/r80/V10/current_realtime_Hop+Onsite_bias_1.0_dt_0.01_rank_80.txt")
dat_HE2 = np.genfromtxt("TT-HEOM/nd2/np10/r100/V10/current_realtime_Hop+Onsite_bias_1.0_dt_0.01_rank_100.txt")
dat_HE3 = np.genfromtxt("TT-HEOM/nd2/np15/r100/V10/current_realtime_Hop+Onsite_bias_1.0_dt_0.01_rank_100.txt")
dat_HE4 = np.genfromtxt("TT-HEOM/nd2/np25/r100/V10/current_realtime_Hop+Onsite_bias_1.0_dt_0.01_rank_100.txt")
dat_HE5 = np.genfromtxt("TT-HEOM/nd2/np25/r150/V10/current_realtime_Hop+Onsite_bias_1.0_dt_0.01_rank_150.txt")

dat_MF = np.genfromtxt("MF-HEOM/V1.0/L3_V1.0_dt0.001_currents.dat")



#plt.plot(dat_ED[:,0],dat_ED[:,1],color="purple",linewidth=7,label=r"ED $(M=6)$")
for ax in [ax0,ax1]:
	ax.plot(dat_MF[:,0],dat_MF[:,1],color="k",linewidth=7,label=r"MF-HEOM")
	ax.plot(dat_HE1[:,0],dat_HE1[:,1],color="orange",linewidth=4,linestyle="dashed",dashes=(9,7),label=r"TT-HEOM [$N_p=10$, $r=80$]")
	ax.plot(dat_HE2[:,0],dat_HE2[:,1],color="purple",linewidth=4,linestyle="dashed",dashes=(9,10),label=r"TT-HEOM [$N_p=10$, $r=100$]")
	ax.plot(dat_HE3[:,0],dat_HE3[:,1],color="blue",linewidth=4,linestyle="dashed",dashes=(10,10),label=r"TT-HEOM [$N_p=15$, $r=100$]")

	ax.plot(dat_HE4[:,0],dat_HE4[:,1],color="red",linewidth=3,linestyle="dashed",dashes=(10,3),label=r"TT-HEOM [$N_p=25$, $r=100$]")
	ax.plot(dat_HE5[:,0],dat_HE5[:,1],color="green",linewidth=3,linestyle="dashed",dashes=(13,6),label=r"TT-HEOM [$N_p=25$, $r=150$]")

#plt.plot(dat_MF[:,0],dat_MF[:,1],color="turquoise",linewidth=4,linestyle="dashed",dashes=(5,10),label="MF-HEOM")

#plt.xlim(0,25)
ax0.set_ylim(-0.075,0.125)
ax0.set_xlabel(r"$t \, [t_0^{-1}]$")
ax0.set_ylabel(r"$I_L\, [t_0]$")
ax1.set_xlim(0,45)
ax1.set_ylim(0.015,0.055)
ax0.legend()
plt.savefig("Comparison_TTHEOM_TDVP_ED_hop+onsite+SOC_V10.pdf")
