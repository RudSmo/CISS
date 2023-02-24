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
fac_t = 1e-3 * 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.055
py = 0.12
sx = 0.27
sy = 0.85
dx = 0.05

fig = plt.figure(figsize=(12,5))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax3 = plt.axes([px+2*(sx+dx)+0.02,py,sx,sy])
axes = [ax1,ax2,ax3]

V = [0,12,20]
Vlat = [L"$V=0\, \mathrm{V}$",L"$V=0.48\, \mathrm{V}$",L"$V=0.8\,\mathrm{V}$"]
w = 0.1 
for (v,vl,ax) in zip(V,Vlat,axes)
	datQ = readdlm("V$(v)_t1_0.004_w0_0.002_nv50.dat")
	ax.plot(datQ[2:end,1]*fac_t,datQ[2:end,2],linewidth=lw,color="k",label="Quantum")
	
	datC = readdlm("../../MFHEOM/eps_-240_w0_2_t0_40_t1_4_avg/V$(v)/phonons_avg_n30.dat")
	ax.plot(datC[2:end,1]*fac_t,datC[2:end,2],linewidth=lw-1,color="red",label="Classical")	
	
	ax.set_xlabel(L"$t\, [\mathrm{ps}]$")
	ax.set_xlim(0,12)
	ax.annotate(vl,(0.6,0.9),xycoords="axes fraction",fontsize=fs)
end

ax1.legend(loc="lower right")
ax1.set_ylabel(L"$\langle X\rangle$",labelpad=-5)
ax1.set_ylim(-0.1,0.5)
ax2.set_ylim(-0.1,50.0)
ax3.set_ylim(-0.2,1.5)


plt.savefig("Peierls_L2_FranssonParams_PhQ.pdf")
