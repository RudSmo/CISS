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

colors = ["k","red","blue"]
size = [3,4,8,12,16]
ms = 12
V = [1,2,3]
Vn = ["0.384","0.44","0.52"]
datname = "Magnetoresistance_time_avg.dat"
for (v,vn,col) in zip(V,Vn,colors)

	dats = ["L4x4/$(datname)",
		  "L3x4/$(datname)",
		  "L2x4/$(datname)",
		  "L2x2/$(datname)", 
		  "L1x3/$(datname)"]


	M_R = []
	s_M_R = []

	for dat in dats[end:-1:1]
		data = readdlm(dat)
		push!(M_R,data[v,2])
		push!(s_M_R,data[v,3])
	end

	
	ax.errorbar(size,M_R,yerr=s_M_R,marker="o",linewidth=2,color=col,label=L"$V=$"*"$(vn)"*L"$\, \mathrm{V}$",markersize=ms,linestyle="")

end

ax.set_xlabel("system size")
ax.set_ylabel(L"$M_R\,[\%]$")
ax.set_xlim(1,17)
ax.legend(loc="lower left")
plt.savefig("Magnetocurrents_sizeVar_P$(P)_VVar.png",dpi=200)
