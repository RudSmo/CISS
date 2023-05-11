using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002 * 4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.7
sx = 0.87
sy = 0.28
dy = 0.03

fig = plt.figure(figsize=(9,9))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])

datMF = readdlm("current_voltage_MFHEOM.dat") 
datLB = readdlm("current_voltage_LB.dat")
datCo = readdlm("conv_voltage.dat") 

axes = [ax1,ax2,ax3]
col = ["blue","green","red"]

Gamma = 0.01
T0 = 0.04
ax1.plot(datLB[2:end,1]/T0,datLB[2:end,2]/Gamma,label="Landauer-Büttiker",linewidth=lw,color=col[3])
ax1.plot(datMF[:,1]/T0,datMF[:,2]/Gamma,label="MFHEOM RK4",linewidth=lw,color="k",linestyle="dashed")
ax1.set_xlim(0,0.35/T0)
ax1.set_xticklabels([])
annx = 0.02
anny = 0.9
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel(L"$I \, [\Gamma]$")
ax1.legend(loc="center right",framealpha=0.8)

ax2.plot(datCo[:,1]/T0,(datCo[:,2] .- 3)*1e7,linewidth=lw,color="k")
ax2.set_xlim(0,0.35/T0)
ax2.set_ylabel(L"$N-\mathrm{tr}\{\rho_\mathrm{steady-state}\} \, [\times 10^{-7}]$")
ax2.set_xticklabels([])
ax2.annotate("(b)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax2.set_ylim(-1,1)

ax3.plot(datCo[:,1]/T0,datCo[:,3],linewidth=lw,color="k")
ax3.set_xlim(0,0.35/T0)
ax3.set_ylabel(L"$E_\mathrm{sys}^\mathrm{steady-state}$")
ax3.annotate("(c)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax3.set_xlabel(L"$V\,[t_0]$")

plt.savefig("L3_uncoupled_open.pdf")

