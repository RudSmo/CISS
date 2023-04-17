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
num_sites = 4
num_laps = 2
N = num_laps * num_sites
V = 0.04 * 0.2 * 48
Vs = ["V48","V55","V65"]
labels = [L"$V=\,0.384\mathrm{V}$",L"$V=\,0.44\mathrm{V}$",L"$V=\,0.52\mathrm{V}$"]

for (v,lab) in zip(Vs,labels)
        dat = readdlm("$(v)/currents_avg_n201.dat")
        ax.plot(dat[:,1]*fac_t*1e-3,(dat[:,2]-dat[:,3])*fac_c,linewidth=lw,label=lab)
end
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
ax.legend()


plt.savefig("L$(num_laps)x$(num_sites)_P0.5_VVar_current.pdf")
