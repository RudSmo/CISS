import numpy as np
import matplotlib.pyplot as plt
import os,fnmatch
from matplotlib import rc

plt.style.use("classic")

ls = 18
lw = 3
rc("font",family="serif",serif="Times")
rc("text",usetex=True)
rc("xtick",labelsize=ls)
rc("ytick",labelsize=ls)
rc("axes",labelsize=ls)
rc("legend",fontsize=ls)


fig = plt.figure(figsize=(10,8))
px = 0.10
py = 0.57
sx = 0.88
sy = 0.41
dy = 0.06

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

elDat = np.genfromtxt("electrons_Vvar.dat")
cuDat = np.genfromtxt("currents_Vvar.dat")

ax1.plot(elDat[:,0],elDat[:,1],linewidth=lw,color="red",label=r"$P_{\uparrow}$")
ax1.plot(elDat[:,0],elDat[:,2],linewidth=lw,color="blue",label=r"$P_{\downarrow}$")
ax1.set_xlim(-0.01,0.51)
ax1.set_ylim(0.0,0.51)
ax1.legend(loc="lower right")
ax1.tick_params(labelbottom=False)
ax1.set_ylabel(r"populations",labelpad=22)

ax2.plot(cuDat[:,0],cuDat[:,1]-cuDat[:,2],linewidth=lw,color="k")
ax2.set_xlim(-0.01,0.51)
ax2.set_ylim(-0.0003,0.01)
ax2.set_ylabel(r"net current $I=I_L-I_R$")
ax2.set_xlabel(r"bias voltage $V$")


plt.savefig("Single_NI_SamParams.pdf")





