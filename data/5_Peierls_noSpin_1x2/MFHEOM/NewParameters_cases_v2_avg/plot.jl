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

plot = "el" # "cur" "ph" "el"

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

params = [("4.0","0.0"),("20.0","0.0"),("2.0","-1.0"),("2.0","-2.0")]
folders = ["i","ii","iii","iv"]
w = 0.1 
colors = ["blue","green","red","orange"]#,"purple","black"]
for (par,fold,col) in zip(params,folders,colors)
	lab = "("*fold*")"
	if plot == "cur"
		cur = readdlm("$(fold)/currents_avg_n30.dat")
		ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=lab)
	elseif plot == "ph"
		cur = readdlm("$(fold)/phonons_avg_n30.dat")
		ax.plot(cur[2:end,1]*fac_t,(w^2/2 .* (cur[2:end,2] .^ 2 .+ cur[2:end,3] .^ 2 ./ w^2)),linewidth=lw,label=lab)
	elseif plot == "el"	
		cur = readdlm("$(fold)/electrons_es_avg_n30.dat")
        	ax.plot(cur[2:end,1]*fac_t,cur[2:end,2],linewidth=lw,label=lab,color=col)
		ax.plot(cur[2:end,1]*fac_t,cur[2:end,3],linewidth=lw,linestyle="dashed",color=col)
	end
end
ax.legend(loc="upper left")

ax.set_xlabel(L"$t\, [\mathrm{fs}]$")

if plot == "cur"
	ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
	plt.savefig("Peierls_L2_NewParams_Cases_Vvar_currents_avg.pdf")
elseif plot == "ph"
	ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{eV}]$")
	ax.set_yscale("log")
	plt.savefig("Peierls_L2_NewParams_Cases_Vvar_phonons_avg.pdf")
elseif plot == "el"
	ax.set_ylabel("Eigenstate population")
	plt.savefig("Peierls_L2_NewParams_Cases_Vvar_ES_avg.pdf")
end
