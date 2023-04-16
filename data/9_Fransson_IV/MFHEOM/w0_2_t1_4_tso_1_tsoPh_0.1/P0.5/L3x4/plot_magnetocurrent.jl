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
num_laps = 3
N = num_laps * num_sites
V = 0.04 * 0.2 * 48
Vs = ["V48","V55","V65"]
labels = [L"$V=\,0.384\mathrm{V}$",L"$V=\,0.44\mathrm{V}$",L"$V=\,0.52\mathrm{V}$"]
for (v,lab) in zip(Vs,labels)
	datP = readdlm("$(v)/currents_avg_n101.dat")
	datM = readdlm("../../P-0.5/L3x4/$(v)/currents_avg_n20.dat")
	curP = datP[:,2]-datP[:,3]
	curM = datM[:,2]-datM[:,3]
	curT = curP + curM
	curD = curP - curM
	ax.plot(datP[:,1]*fac_t*1e-3,100 * curD ./ curT,linewidth=lw,label=lab)
end


ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$M_R\,[\%]$")
ax.legend()


plt.savefig("L$(num_laps)x$(num_sites)_VVar_magnetocurrent.pdf")
