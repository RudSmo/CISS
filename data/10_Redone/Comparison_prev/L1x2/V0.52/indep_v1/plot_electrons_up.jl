using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+sx+dx,py,0.02,sy])
num_sites = 2
num_laps = 1
N = num_laps * num_sites
vcur = 65
V = 0.04 * 0.2 * vcur
dat = readdlm("electrons_avg_n41.dat")
cp = plt.cm.jet
cmap = plt.cm.ScalarMappable(cmap=cp)

vmn = minimum(dat[:,2:2:end])
vmx = maximum(dat[:,2:2:end])

ax.imshow(dat[:,2:2:end],aspect="auto",vmin=vmn,vmax=vmx,interpolation="none",extent=[1,N,0,dat[end,1]*fac_t],cmap="jet")
ax.set_xlabel(L"site $i$")
ax.set_ylabel(L"$t \, [\mathrm{fs}]$",labelpad=-0.5)
cbar = fig.colorbar(cmap, ax=ax, cax = cax, aspect=90,pad=.001)
cbar.set_label(L"$n_{i\uparrow}$")
cbar.ax.tick_params(labelsize=fs)


plt.savefig("L$(num_laps)x$(num_sites)_V$(V)_electrons_up.pdf")
