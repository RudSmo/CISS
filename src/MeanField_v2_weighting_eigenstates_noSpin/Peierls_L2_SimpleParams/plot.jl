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

V = ["3.0"]#,"4.0","7.0"]
Vlat = [L"$3.0$"]#,L"$4.0$",L"$7.0$"]
w = 0.1 
colors = ["blue","green","red"]
for (v,vl,col) in zip(V,Vlat,colors)
#	cur = readdlm("L2_gPh0.0_tPh0.2_V$(v)_w00.1_dt0.0005_1_currents.dat")
#	ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=L"$V=\,$"*vl*" V")

#	cur = readdlm("L2_gPh0.0_tPh0.2_V$(v)_w00.1_dt0.0005_1_phonons.dat")
#	ax.plot(cur[2:end,1]*fac_t,(w^2/2 .* (cur[2:end,2] .^ 2 .+ cur[2:end,3] .^ 2 ./ w^2)),linewidth=lw,label=L"$V=\,$"*vl*" V")
	
	cur = readdlm("L2_gPh0.0_tPh0.2_V$(v)_w00.1_dt0.0005_1_electrons_es.dat")
        ax.plot(cur[2:end,1]*fac_t,cur[2:end,2],linewidth=lw,label=L"$V=\,$"*vl*" V",color=col)
	ax.plot(cur[2:end,1]*fac_t,cur[2:end,3],linewidth=lw,linestyle="dashed",color=col)
end
ax.legend()
ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
#ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")

#ax.set_xlim(0,3000)
#ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{eV}]$")

ax.set_ylabel("Eigenstate population")

#plt.savefig("Peierls_L2_NewParams_Vvar_phonons.pdf")
#plt.savefig("Peierls_L2_NewParams_Vvar.pdf")
plt.savefig("Peierls_L2_NewParams_V3.0_ES.pdf")
