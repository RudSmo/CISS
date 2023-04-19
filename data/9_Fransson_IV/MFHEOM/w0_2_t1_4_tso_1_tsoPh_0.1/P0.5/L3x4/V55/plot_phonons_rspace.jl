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
sx = 0.78
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(9,7))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+sx+dx,py,0.01,sy])
num_sites = 4
num_laps = 3
tPh = 0.1 * 0.04 
ws = 0.002
fac_X = 1.9732*1e-7 * 10^10 * sqrt(1e-9)
N = num_laps * num_sites
vcur = 55
V = 0.04 * 0.2 * vcur
P = 0.5
dat = readdlm("phonons_avg_n51.dat")
cp = plt.cm.inferno
cmap = plt.cm.ScalarMappable(cmap=cp)
sp = dat[:,2:2:end-2] * fac_X
vmn = minimum(sp)
vmx = maximum(sp)


im = ax.imshow(sp,aspect="auto",vmin=vmn,vmax=vmx,interpolation="none",extent=[1,N,dat[end,1]*fac_t*1e-3,0],cmap="inferno")
ax.set_xlabel(L"site $i$")
ax.set_ylabel(L"$t \, [\mathrm{ps}]$",labelpad=7)
cbar = fig.colorbar(im, ax=ax, cax = cax, aspect=90,pad=.001)
cbar.set_label("inter-site bond displacement "*L"$X \, [\frac{\r{A}}{\sqrt{\mathrm{GeV}}}]$")
cbar.ax.tick_params(labelsize=fs)


plt.savefig("L$(num_laps)x$(num_sites)_V$(V)_P$(P)_phonons_rspace.pdf")
