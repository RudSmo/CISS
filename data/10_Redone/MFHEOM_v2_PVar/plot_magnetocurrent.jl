using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.11
py = 0.12
sx = 0.87
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(9,4.5))
ax = plt.axes([px,py,sx,sy])
P = 0.5

colors = ["k","red","blue"]
symbols = ["x","o","d"]

#syssizes = [["L1x3","L2x3","L3x3"],
#	    ["L2x4","L3x4"],
# 	    ["L2x2","L2x5","L2x6"]]
#
#sizes = [[3,6,9],
#	 [8,12],
#	 [4,10,12]]

syssizes = [["L1x3","L2x2","L2x3","L2x4","L3x3","L2x5","L3x4"]]
sizes = [[3,4,6,8,9,10,12]]

ms = 9
V = [1,2,3]
Vn = ["0.4","0.48","0.6"]

datname = "Magnetoresistance_time_avg_stderr.dat"


for (v,vn,col) in zip(V,Vn,colors)
	for (size,sys,symb) in zip(sizes,syssizes,symbols)



		M_R = []
		s_M_R = []

		for s in sys
			
			data = readdlm("$(s)/$(datname)")
			push!(M_R,abs(data[v,2]))
			push!(s_M_R,data[v,3])
		end

		if size == sizes[1]
			ax.errorbar(size,M_R,yerr=s_M_R,marker=symb,linewidth=2,color=col,label=L"$V=\,$"*"$(vn)"*L"$\, \mathrm{V}$",markersize=ms)
		else
			ax.errorbar(size,M_R,yerr=s_M_R,marker=symb,linewidth=2,color=col,markersize=ms)
		end
	end
end

ax.set_xlabel("system size "*L"$N$")
ax.set_ylabel(L"$M_R\,[\%]$")
ax.set_xlim(2,13)
ax.legend(loc="upper left")
#plt.savefig("Magnetocurrents_sizeVar_VVar.png",dpi=200)
plt.savefig("Magnetocurrents_sizeVar_VVar_stderr.pdf")
