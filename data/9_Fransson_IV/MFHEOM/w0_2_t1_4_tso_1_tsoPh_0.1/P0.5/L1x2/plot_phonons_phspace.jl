using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)
ws = 0.002
fac_X = 1.9732*1e-7 * 1e10 * ws
fac_P = 5.3442883*1e-28 * 1e-12 / (sqrt(ws) * 6.022*1e26 * 1e-10)
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.12
sx = 0.88
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
num_sites = 2
num_laps = 1
N = num_laps * num_sites
vcur = 48
V = 0.04 * 0.2 * vcur
P = 0.5
dat = readdlm("phonons_avg_n41.dat")
for i in 1:2:Int(N)
        ax.plot(dat[:,i+1]*fac_X,dat[:,i+2]*fac_P,linewidth=lw,label="Bond $((i+1)/2)")	
end
ax.set_xlabel(L"$\omega_0 X\, [\r{A}]$")
ax.set_ylabel(L"$P\, \left[\frac{\mathrm{amu}\cdot \r{A}}{\mathrm{ps}}\right]$")
ax.legend()

plt.savefig("L$(num_laps)x$(num_sites)_P$(P)_V$(V)_phonons.pdf")
