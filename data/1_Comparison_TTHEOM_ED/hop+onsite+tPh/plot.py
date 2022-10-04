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


fig = plt.figure(figsize=(10,8))

dat_ED = np.genfromtxt("ED/office/nvib6/pops_nvib6_vib.dat")
dat_HE1 = np.genfromtxt("HEOM/r100/n_site_realtime_bias_0.0_dt_0.001_rank_100.txt")
dat_HE2 = np.genfromtxt("HEOM/r150/n_site_realtime_bias_0.0_dt_0.001_rank_150.txt")
dat_HEn1 = np.genfromtxt("HEOM/r100/norm_realtime_bias_0.0_dt_0.001_rank_100.txt")
dat_HEn2 = np.genfromtxt("HEOM/r150/norm_realtime_bias_0.0_dt_0.001_rank_150.txt")
#dat_MF = np.genfromtxt("")
dat_MP = np.genfromtxt("tdvp/N3_t01.00_tso0.000_tPh0.100_dt0.00500.dat")

plt.plot(dat_ED[:,0],dat_ED[:,1],color="k",linewidth=7,label=r"ED $(M=6)$")
plt.plot(dat_MP[:,0],dat_MP[:,1],color="red",linewidth=4,linestyle="dashed",dashes=(7,5),label=r"TDVP $(M=5)$")
plt.plot([0],[0],color="white",label="TT-HEOM:")
plt.plot(dat_HE1[:,0],dat_HE1[:,1],color="gray",linewidth=4,linestyle="dashed",dashes=(6,7),label=r"$(M=5), r=100$")
plt.plot(dat_HEn1[:,0],dat_HEn1[:,1],color="gray",linewidth=4,label=r"$(M=5), r=100$ $\vert \vert\rho\rangle\vert$")
plt.plot(dat_HE2[:,0],dat_HE2[:,1],color="blue",linewidth=4,linestyle="dashed",dashes=(6,9),label=r"$(M=5), r=150$")
plt.plot(dat_HEn2[:,0],dat_HEn2[:,1],color="blue",linewidth=4,label=r"$(M=5), r=150$ $\vert \vert\rho\rangle\vert$")
#plt.plot(dat_MF[:,0],dat_MF[:,1],color="turquoise",linewidth=4,linestyle="dashed",dashes=(5,10),label="MF-HEOM")

plt.xlim(0,25)
plt.ylim(-0.1,1.25)
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$n_{1\uparrow}$, $\vert \vert\rho\rangle\vert $")

plt.legend()
plt.savefig("Comparison_TTHEOM_TDVP_ED_hop+onsite+tph.pdf")
