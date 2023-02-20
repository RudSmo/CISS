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

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

params = [("4.0","0.0"),("20.0","0.0"),("2.0","-1.0"),("2.0","-2.0")]#,("2.0","1.0"),("2.0","2.0")]
labels = ["(i)","(ii)","(iii)","(iv)"]#,"(v)","(vi)"]
w = 0.1 
colors = ["blue","green","red","orange"]#,"purple","black"]
for (par,lab,col) in zip(params,labels,colors)
	filn = "L2_gPh0.0_tPh0.2_V$(par[1])_w00.1_E$(par[2])_dt0.0005_1_"
#	cur = readdlm("$(filn)currents.dat")
#	ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=lab)

	cur = readdlm("$(filn)phonons.dat")
	ax.plot(cur[2:end,1]*fac_t,(w^2/2 .* (cur[2:end,2] .^ 2 .+ cur[2:end,3] .^ 2 ./ w^2)),linewidth=lw,label=lab)
	
#	cur = readdlm("$(filn)electrons_es.dat")
#        ax.plot(cur[2:end,1]*fac_t,cur[2:end,2],linewidth=lw,label=lab,color=col)
#	ax.plot(cur[2:end,1]*fac_t,cur[2:end,3],linewidth=lw,linestyle="dashed",color=col)
end
ax.legend(loc="upper left")

ax.set_xlabel(L"$t\, [\mathrm{fs}]$")

#ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
#plt.savefig("Peierls_L2_NewParams_Cases_Vvar_currents.pdf")

ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{eV}]$")
ax.set_yscale("log")
plt.savefig("Peierls_L2_NewParams_Cases_Vvar_phonons.pdf")

#ax.set_ylabel("Eigenstate population")
#plt.savefig("Peierls_L2_NewParams_Cases_Vvar_ES.pdf")
