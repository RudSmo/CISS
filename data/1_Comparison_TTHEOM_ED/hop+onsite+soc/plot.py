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

dat_ED = np.genfromtxt("data_ED_hop+onsite+soc.dat")
dat_HE = np.genfromtxt("data_HEOM_hop+onsite+soc/n_site_r100.dat")
dat_MF = np.genfromtxt("data_MFHEOM/L3_Hop+Onsite+SOC_dt0.001_1_electrons.dat")
dat_MP = np.genfromtxt("data_MPO_hop+onsite+soc.dat")

plt.plot(dat_HE[:,0],dat_HE[:,1],color="k",linewidth=7,label="TT-HEOM")
plt.plot(dat_MF[:,0],dat_MF[:,1],color="red",linewidth=4,linestyle="dashed",dashes=(7,5),label="MF-HEOM")
plt.plot(dat_MP[:,0],dat_MP[:,1],color="gray",linewidth=5,linestyle="dashed",dashes=(6,7),label="MPO-DM")
plt.plot(dat_ED[:,0],dat_ED[:,1],color="turquoise",linewidth=4,linestyle="dashed",dashes=(5,10),label="ED")

plt.xlim(0,25)
plt.ylim(-0.1,1.1)
plt.xlabel(r"$t \, [t_0^{-1}]$")
plt.ylabel(r"$n_{1\uparrow}$")

plt.legend()
plt.savefig("Comparison_TTHEOM_MFHEOM_MPO_ED_hop+onsite+soc.pdf")
