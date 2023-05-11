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
	dat = readdlm("$(size)/$(size)_currents_time_P0.5_avg.dat")
	dat_std = readdlm("$(size)/$(size)_currents_time_P0.5_std.dat")
	volts = dat[:,1]
	M_R = (dat[:,2]-dat[:,3])*fac_c
	s_M_R = sqrt.(dat_std[:,2].^2 + dat_std[:,3].^2) * fac_c
	ax.errorbar(volts,M_R,yerr=s_M_R,marker="o",linewidth=2,label=L"$L=\,$"*"$(n)"*L"$\, \times 4$",markersize=ms)

end

ax.set_xlabel("voltage "*L"$V\,[\mathrm{V}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
#ax.set_xlim(1,17)
ax.set_ylim(0.0,1.5)
ax.legend(loc="lower left")
plt.savefig("Currents_1Mode_sizeVar_VVar.png",dpi=200)
#plt.savefig("Currents_1Mode_sizeVar_VVar.pdf",dpi=200)
