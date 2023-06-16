using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.05
py = 0.12
sx = 0.28
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(13,7))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
V = 60
cax = plt.axes([px+3*(sx+dx),py,0.01,sy])

dat1 = readdlm("L1x3/V60/electrons_avg_n61_P.dat")
dat2 = readdlm("L2x3/V60/electrons_avg_n71_P.dat")
dat3 = readdlm("L3x3/V60/electrons_avg_n19_P.dat")

cp = plt.cm.inferno
cmap = plt.cm.ScalarMappable(cmap=cp)
sp1 = dat1[2:end,2:2:end] .- dat1[2:end,3:2:end]
sp2 = dat2[2:end,2:2:end] .- dat2[2:end,3:2:end]
sp3 = dat3[2:end,2:2:end] .- dat3[2:end,3:2:end]
vmn = minimum(sp2)
vmx = maximum(sp2)-0.05

N = 3
ax1.imshow(sp1,aspect="auto",vmin=vmn,vmax=vmx,cmap="inferno",extent=[0.5,N+0.5,dat1[end,1]*fac_t*1e-3,0])
N = 6
ax2.imshow(sp2,aspect="auto",vmin=vmn,vmax=vmx,extent=[0.5,N+0.5,dat2[end,1]*fac_t*1e-3,0],cmap="inferno")
N = 9
im = ax3.imshow(sp3,aspect="auto",vmin=vmn,vmax=vmx,extent=[0.5,N+0.5,dat3[end,1]*fac_t*1e-3,0],cmap="inferno")

ax1.set_xlabel(L"site $i$")
ax1.annotate(L"$L=1\times 3$",(0.7,0.05),xycoords="axes fraction",fontsize=fs,color="white")
ax1.set_ylabel(L"$t \, [\mathrm{ps}]$",labelpad=7)

ax2.set_xlabel(L"site $i$")
ax2.annotate(L"$L=2\times 3$",(0.7,0.05),xycoords="axes fraction",fontsize=fs,color="white")

ax3.set_xlabel(L"site $i$")
ax3.annotate(L"$L=3\times 3$",(0.7,0.05),xycoords="axes fraction",fontsize=fs,color="white")

ax1.set_ylim(85,0)
ax2.set_yticklabels([])
ax2.set_ylim(85,0)
ax3.set_yticklabels([])
ax3.set_ylim(85,0)

cbar = fig.colorbar(im, ax=ax3, cax = cax, aspect=90,pad=.001)
cbar.set_label(L"$s_{i}^z$")
cbar.ax.tick_params(labelsize=fs)

plt.savefig("Sz_space_dyn_sysVar_x3.pdf")
