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

V = ["0.0","0.4","0.48","0.56","1.0"]
Vlat = [L"$0.0$",L"$0.4$",L"$0.48$",L"$0.56$",L"$1.0$"]
w = 2.0 
for (v,vl) in zip(V,Vlat)
	cur = readdlm("L2_gPh0.0_tPh-0.004_V$(v)_w00.002_E-0.24_dt0.05_1_currents.dat")
	ax.plot(cur[2:end,1]*fac_t,(cur[2:end,2]-cur[2:end,3])*fac_c,linewidth=lw,label=L"$V=\,$"*vl*" V")

#	cur = readdlm("L2_gPh0.0_tPh-0.004_V$(v)_w00.002_E-0.24_dt0.05_1_phonons.dat")
#	ax.plot(cur[2:end,1]*fac_t, 1e-3*(w^2/2 .* cur[2:end,2] .^ 2 .+ 0.5 * cur[2:end,3] .^ 2),linewidth=lw,label=L"$V=\,$"*vl*" V")
end
ax.legend()

ax.set_xlim(0,12)

ax.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
plt.savefig("Peierls_L2_FranssonParams_Vvar_currents.pdf")

#ax.set_ylabel(L"$E_\mathrm{ph}\,[\mathrm{eV}]$")
#plt.savefig("Peierls_L2_FranssonParams_Vvar_phonons.pdf")

