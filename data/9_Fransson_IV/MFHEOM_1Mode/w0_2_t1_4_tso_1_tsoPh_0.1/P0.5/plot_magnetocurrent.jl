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

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.09
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(9,7))
ax = plt.axes([px,py,sx,sy])
P = 0.5

sizes = ["L2x4","L3x4","L4x4","L5x4"]
ns = [2,3,4,5]

ms = 12

for (size,n) in zip(sizes,ns)
	dat = readdlm("$(size)/$(size)_magnetoresistance_avg_std.dat")
	volts = dat[:,1]
	M_R = dat[:,2]
	s_M_R = dat[:,3]
	ax.errorbar(volts,M_R,yerr=s_M_R,marker="o",linewidth=2,label=L"$L=\,$"*"$(n)"*L"$\, \times 4$",markersize=ms)

end

ax.set_xlabel("voltage "*L"$V\,[\mathrm{V}]$")
ax.set_ylabel(L"$M_R\,[\%]$")
#ax.set_xlim(1,17)
ax.set_ylim(-15,15)
ax.legend(loc="lower left")
plt.savefig("Magnetocurrents_1Mode_sizeVar_VVar.png",dpi=200)
#plt.savefig("Magnetocurrents_1Mode_sizeVar_VVar.pdf",dpi=200)
