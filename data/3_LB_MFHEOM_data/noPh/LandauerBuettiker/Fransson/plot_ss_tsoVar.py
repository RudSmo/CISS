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
ax1 = plt.axes([0.12,0.147,0.85,0.8])
#ax1 = plt.axes([0.1,0.5,0.86,0.45])
#ax2 = plt.axes([0.1,0.12,0.86,0.45])


dat1P = np.genfromtxt("7x8/fS8_nL7_tso0.001_P0.5_currents_v3.dat")
dat1M = np.genfromtxt("7x8/fS8_nL7_tso0.001_P-0.5_currents_v3.dat")

dat2P = np.genfromtxt("7x8/fS8_nL7_tso0.05_P0.5_currents_v3.dat")
dat2M = np.genfromtxt("7x8/fS8_nL7_tso0.05_P-0.5_currents_v3.dat")
dat3P = np.genfromtxt("7x8/fS8_nL7_tso0.01_P0.5_currents_v3.dat")
dat3M = np.genfromtxt("7x8/fS8_nL7_tso0.01_P-0.5_currents_v3.dat")
Gamma = 0.25

def SpinSelectivity(dataP,dataM):
	SP = []
	V = []
	for j in range(len(dataP[:,0])):	
		if dataP[j,1] > 1e-5 and dataM[j,1] > 1e-5:
			sp = 100 * (dataP[j,1]-dataM[j,1]) / (dataP[j,1]+dataM[j,1])
			V.append(dataP[j,0])
			SP.append(sp)
#		else:
#			sp = 0.0
	return V,np.array(SP)

V1,SP1 = SpinSelectivity(dat1P,dat1M)
V2,SP2 = SpinSelectivity(dat2P,dat2M)
V3,SP3 = SpinSelectivity(dat3P,dat3M)

ax1.plot(V1,SP1*1e5,linewidth=6,label=r"$L=7\times 8, \lambda_0 = 1\,\mathrm{meV}$")
ax1.plot(V3,SP3*1e5,linewidth=6,label=r"$\lambda_0 = 10\,\mathrm{meV}$")
ax1.plot(V2,SP2*1e5,linewidth=6,label=r"$\lambda_0 = 50\,\mathrm{meV}$")
ax1.legend()


#ax1.set_xlim(-0.1,20)
ax1.set_ylabel(r"$S\, [\times 10^{-5} \%]$")
ax1.set_xlabel(r"$V \, [\mathrm{V}]$")

plt.savefig("LB_tsoVar_SP.pdf")
