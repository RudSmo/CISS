using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


fig = plt.figure(figsize = (10,4.5))
px = 0.08
py = 0.12
sx = 0.89
dx = 0.1
sy = 0.84
ax = plt.axes([px,py,sx,sy])


d1 = readdlm("L2_tPh0.2_SOC0.0_SOCPh0.0_V0.0_w00.1_P0.0_dt0.001_1_electrons.dat")
d2 = readdlm("gL0.25_gR0.25_w0_0.1_tph0.2_electrons.dat")


ax.plot(d1[2:end,1],d1[2:end,2],color="k",linewidth=3,label="MH-HEOM")
ax.plot(d2[2:end,1],d2[2:end,2],color="red",linewidth=2,linestyle="dashed",label=L"Liouville-von Neumann under $\mathbf{H}^\mathrm{eff}$")
ax.set_xlabel(L"$t \, [t_0^{-1}]$",labelpad=-3)
ax.set_ylabel(L"$n$")
ax.legend(loc="upper right")
#ax.set_xlim(0,50)

plt.savefig("Comparison_NH_HEOM_withPh.pdf")



