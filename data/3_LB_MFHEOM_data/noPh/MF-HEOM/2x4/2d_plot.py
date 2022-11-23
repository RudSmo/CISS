import numpy as np
import matplotlib.pyplot as plt
import os,fnmatch
from matplotlib import rc
from matplotlib import cm
import matplotlib.colors as colors

plt.style.use("classic")

ls = 23
lw = 3

rc("font",family="serif",serif="Times")
rc("text",usetex=True)
rc("xtick",labelsize=ls)
rc("ytick",labelsize=ls)
rc("axes",labelsize=ls)
rc("legend",fontsize=ls-5)


fig = plt.figure(figsize=(10,11))
px = 0.08
py = 0.69
sx = 0.78
sy = 0.29
dy = 0.02
dx = 0.01

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
axc = plt.axes([px+sx+dx,py-2*(sy+dy),0.02,3*sy+2*dy])

Ps = [0.5,0.0,-0.5]
axes = [ax1,ax2,ax3]
vmin = -0.25
vmax = 0.25
cp = plt.cm.RdBu_r
norm = colors.SymLogNorm(linthresh=0.002, linscale=1.0,vmin=vmin,vmax=vmax, base=10)
cmap = cm.ScalarMappable(norm=norm, cmap=cp)


for P,ax in zip(Ps,axes):
	Pf = f"P{P}/"
	datE = np.genfromtxt(Pf+"Electrons_steadyState.dat")
	im = ax.imshow(datE[::-1,1::2]-datE[::-1,2::2],aspect="auto",extent=[1,8,0,datE[-1,0]],norm=norm,cmap = "RdBu_r")
	ax.annotate(rf"$P={P}$",(0.8,0.1),xycoords="axes fraction",fontsize=ls)
		
	ax.set_ylim(0,15)
	if ax != ax3:
		ax.set_xticklabels([])
	ax.set_ylabel(r"$V\, [t_0]$")
	
cbar = fig.colorbar(cmap, ax=ax1, cax = axc, aspect=90,pad=.001)
cbar.set_label(r"$n_{i\uparrow}-n_{i\downarrow}$")
cbar.ax.tick_params(labelsize=17)

ax3.set_xlabel(r"sites")
plt.savefig("MF-HEOM_2x4_Sz_steady_PVar.pdf")
