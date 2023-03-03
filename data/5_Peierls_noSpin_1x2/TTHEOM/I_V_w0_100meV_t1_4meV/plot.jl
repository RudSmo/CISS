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

V = [0,10,12,14]
Vlat = [L"$V=0\,\mathrm{V}$", 
	L"$V=0.4\,\mathrm{V}$",
	L"$V=0.48\,\mathrm{V}$",
	L"$V=0.56\,\mathrm{V}$",
	L"$V=0.8\,\mathrm{V}$"]
w0 = 0.1
t1 = 0.004 
colors = ["blue","green","red","orange","gray","magenta"]



for (v,vl,col) in zip(V,Vlat,colors)
	fname = "V$(v)/"
	curMF = readdlm("../../MFHEOM/eps_-240_w0_100_t1_4/$(fname)currents_avg_n30.dat")
	ax.plot(curMF[2:end,1]*fac_t,(curMF[2:end,2]-curMF[2:end,3])*fac_c*2,linewidth=lw,label=vl,color=col)
	
	ve = v * 0.04
	curTT = readdlm("$(fname)current_realtime_Hop+Onsite_bias_$(ve)_dt_0.25_rank_100.txt")
	ax.plot(curTT[2:end,1]*fac_t,(curTT[2:end,2]-curTT[2:end,3])*fac_c,linewidth=lw-1,color="k",linestyle="dashed")
end

ax.legend(loc="center right")
ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
ax.set_xlim(0,750)
ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
plt.savefig("Peierls_L2_w0_$(w0)_t1_$(t1)_Vvar_comparison_currents.pdf")
