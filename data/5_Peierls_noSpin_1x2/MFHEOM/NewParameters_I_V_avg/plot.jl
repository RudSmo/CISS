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

obs = "ph" # "cur"  "ph", "el"

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

V = ["2.0","4.0","4.2","4.4","4.6","5.0"]
Vlat = [L"$2.0$",L"$4.0$",L"$4.2$",L"$4.4$",L"$4.6$",L"$5.0$"]
w = 0.1 
colors = ["blue","green","red","orange","gray","magenta"]



for (v,vl,col) in zip(V,Vlat,colors)
	Vfold = Int(parse(Float64,v)*10)
	fname = "V$(Vfold)/"
	if obs == "cur"
		cur = readdlm("$(fname)currents_avg_n30.dat")
		ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=L"$V=\,$"*vl*" V")

	elseif obs == "ph"
		cur = readdlm("$(fname)phonons_avg_n30.dat")
		ax.plot(cur[2:end,1]*fac_t,w^2/2 .* cur[2:end,2] .^ 2 .+ cur[2:end,3] .^ 2 ./ 2,linewidth=lw,label=L"$V=\,$"*vl*" V")
	elseif obs == "el"	
		cur = readdlm("$(fname)electrons_es_avg_n30.dat")
	        ax.plot(cur[2:end,1]*fac_t,cur[2:end,2],linewidth=lw,label=L"$V=\,$"*vl*" V",color=col)
		ax.plot(cur[2:end,1]*fac_t,cur[2:end,3],linewidth=lw,linestyle="dashed",color=col)
	end
end
ax.legend(loc="center right")
ax.set_xlabel(L"$t\, [\mathrm{fs}]$")

if obs == "cur"
	ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
	plt.savefig("Peierls_L2_NewParams_Vvar_avg_currents.pdf")
elseif obs == "ph"
	ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{eV}]$")
	plt.savefig("Peierls_L2_NewParams_Vvar_avg_phonons.pdf")
elseif obs == "el"
	ax.set_ylabel("Eigenstate population")
	plt.savefig("Peierls_L2_NewParams_Vvar_avg_ES.pdf")
end
