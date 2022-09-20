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

dat_ED1 = np.genfromtxt("ED/M5/pops_vib.dat")
dat_ED2 = np.genfromtxt("ED/M10/pops_vib.dat")
dat_HE1 = np.genfromtxt("heom/M5_r100/n_site_realtime_bias_0.0_dt_0.001_rank_100.txt")
dat_HE1n = np.genfromtxt("heom/M5_r100/norm_realtime_bias_0.0_dt_0.001_rank_100.txt")
dat_HE2 = np.genfromtxt("heom/M10_r200/n_site_realtime_bias_0.0_dt_0.001_rank_200.txt")
dat_HE2n = np.genfromtxt("heom/M10_r200/norm_realtime_bias_0.0_dt_0.001_rank_200.txt")

plt.plot(dat_ED1[:,0],dat_ED1[:,1],color="k",linewidth=5,label=r"ED $(M=5)$")
plt.plot(dat_ED2[:,0],dat_ED2[:,1],color="red",linewidth=2,linestyle="dashed",label=r"ED $(M=10)$")

#plt.plot(dat_MP[:,0],dat_MP[:,1],color="red",linewidth=4,linestyle="dashed",dashes=(7,5),label=r"TDVP $(M=5)$")
plt.plot([0],[0],color="white",label="TT-HEOM:")
plt.plot(dat_HE1[:,0],dat_HE1[:,1],color="gray",linewidth=3,label=r"$M=5, r=100$")
#plt.plot(dat_HEn1[:,0],dat_HEn1[:,1],color="gray",linewidth=4,label=r"$(M=5), r=100$ $\vert \vert\rho\rangle\vert$")
plt.plot(dat_HE2[:,0],dat_HE2[:,1],color="blue",linewidth=1.5,label=r"$M=10, r=200$")
#plt.plot(dat_HEn2[:,0],dat_HEn2[:,1],color="orange",linewidth=4,linestyle="dashed",dashes=(6,9),label=r"$(M=5), r=200$ $\vert \vert\rho\rangle\vert$")
#plt.plot(dat_MF[:,0],dat_MF[:,1],color="turquoise",linewidth=4,linestyle="dashed",dashes=(5,10),label="MF-HEOM")

plt.xlim(0,200)
plt.ylim(-0.1,1.25)
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$n_{1\uparrow}$, $\vert \vert\rho\rangle\vert $")

plt.legend()
plt.savefig("Comparison_L2_TTHEOM_ED_hop+onsite+soc+tph+socph.pdf")
