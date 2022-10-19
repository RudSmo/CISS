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


fig = plt.figure(figsize=(10,9))

dat_ED = np.genfromtxt("ED/office/nvib6/pops_nvib6_vib.dat")
dat_HE1 = np.genfromtxt("HEOM/Leads_Phonons_spinUncomb/r25/n_site_realtime_Hop+Onsite+tPh_bias_0.0_dt_0.01_rank_25.txt")
dat_HE2 = np.genfromtxt("HEOM/Leads_Phonons_spinUncomb/r80/n_site_realtime_Hop+Onsite+tPh_bias_0.0_dt_0.01_rank_80.txt")
dat_HE3 = np.genfromtxt("HEOM/Leads_Phonons_spinUncomb/r120/n_site_realtime_Hop+Onsite+tPh_bias_0.0_dt_0.01_rank_120.txt")
dat_MP = np.genfromtxt("tdvp/N3_t01.00_tso0.000_tPh0.100_dt0.00500.dat")

plt.plot(dat_ED[:,0],dat_ED[:,1],color="purple",linewidth=7,label=r"ED $(M=6)$")
plt.plot(dat_MP[:,0],dat_MP[:,1],color="orange",linewidth=5,linestyle="dashed",dashes=(7,5),label=r"GSE-TDVP (iTensor) $[M=5, r_\mathrm{max}=47]$")
plt.plot(dat_HE1[:,0],dat_HE1[:,1],color="k",linewidth=5,linestyle="dashed",dashes=(6,7),label=r"TT-HEOM, 1TDVP $[M=5, r=25]$")
plt.plot(dat_HE2[:,0],dat_HE2[:,1],color="red",linewidth=5,linestyle="dashed",dashes=(5,10),label=r"TT-HEOM, 1TDVP $[M=5, r=80]$")
plt.plot(dat_HE3[:,0],dat_HE3[:,1],color="blue",linewidth=5,linestyle="dashed",dashes=(6,7),label=r"TT-HEOM, 1TDVP $[M=5, r=120]$")
#plt.plot(dat_MF[:,0],dat_MF[:,1],color="turquoise",linewidth=4,linestyle="dashed",dashes=(5,10),label="MF-HEOM")

plt.xlim(0,25)
plt.ylim(-0.1,1.25)
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$n_{1\uparrow}$")

plt.legend()
plt.savefig("Comparison_TTHEOM_TDVP_ED_hop+onsite+tph.png")
