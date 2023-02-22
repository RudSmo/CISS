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
fac_t = 1e-3/(4.135668)

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

Vf = [0,10,12,14,20]
V = ["0.0","0.4","0.48","0.56","1.0"]
Vlat = [L"$0.0$",L"$0.4$",L"$0.48$",L"$0.56$",L"$1.0$"]
w = 0.002
for (v,vf,vl) in zip(V,Vf,Vlat)
#	cur = readdlm("V$(vf)/currents_avg_n30.dat")
#	ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=L"$V=\,$"*vl*" V")

	cur = readdlm("V$(vf)/phonons_avg_n30.dat")
	ax.plot(cur[2:end,1]*fac_t, 1e3*(w^2/2 .* cur[2:end,2] .^ 2 .+ 0.5 * cur[2:end,3] .^ 2),linewidth=lw,label=L"$V=\,$"*vl*" V")
end
ax.legend()

ax.set_xlim(0,12)

ax.set_xlabel(L"$t\, [\mathrm{ps}]$")

#ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
#plt.savefig("Peierls_L2_FranssonParams_Vvar_avg_currents.pdf")

ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{meV}]$")
plt.savefig("Peierls_L2_FranssonParams_Vvar_avg_phonons.pdf")

