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
lw = 3
fac_t = 1/(4.135668*1e3) 
fac_c = 1e6 * 2.8494561 

px = 0.08
py = 0.096
sx = 0.89
sy = 0.88
dy = 0.01

fig = plt.figure(figsize=(10,6))
ax = plt.axes([px,py,sx,sy])
ax2 = plt.axes([0.5,0.15,0.4,0.2])


dat1 = readdlm("L8_tPh0.0_SOC0.001_SOCPh0.0_V0.0_w00.002_P0.5_X0_-1.0_dt0.05_1_currents.dat")
dat2 = readdlm("../Comp_SOC/V0.0/L8_tPh0.0_SOC0.001_SOCPh0.0_V0.0_w00.002_P0.0_X0_-1.0_dt0.05_1_currents.dat")
dat3 = readdlm("../Comp_SOC/V0.0/L8_tPh0.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.5_X0_-1.0_dt0.05_1_currents.dat")
dat4 = readdlm("../Comp_SOC/V0.0/L8_tPh0.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_X0_-1.0_dt0.05_1_currents.dat")

ax.plot(dat1[2:end,1]*fac_t,(dat1[2:end,2]-dat1[2:end,3])*fac_c,color="k",linewidth=lw+2,label=L"$P=0.5, \, \lambda_0=0.001\,\mathrm{eV}$")
ax.plot(dat2[2:end,1]*fac_t,(dat2[2:end,2]-dat2[2:end,3])*fac_c,color="red",linewidth=lw+2,label=L"$P=0.0,\, \lambda_0=0.001\,\mathrm{eV}$")
ax.plot(dat3[2:end,1]*fac_t,(dat3[2:end,2]-dat3[2:end,3])*fac_c,color="green",linewidth=lw,label=L"$P=0.5,\, \lambda_0=0.0\,\mathrm{eV}$")
ax.plot(dat4[2:end,1]*fac_t,(dat4[2:end,2]-dat4[2:end,3])*fac_c,color="orange",linestyle="dashed",linewidth=lw,label=L"$P=0.0,\, \lambda_0=0.0\,\mathrm{eV}$")
ax.legend(loc="best")



ax.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax.set_ylabel(L"$I_L - I_R \, [\mathrm{nA}]$")
ax.set_ylim(-20,20)
ax.set_xlim(0,1.8)

ax2.plot(dat1[2:end,1]*fac_t,(dat1[2:end,2]-dat1[2:end,3])*fac_c,color="k",linewidth=lw+2,label=L"$P=0.5, \, \lambda_0=0.001\,\mathrm{eV}$")
ax2.plot(dat2[2:end,1]*fac_t,(dat2[2:end,2]-dat2[2:end,3])*fac_c,color="red",linewidth=lw+2,label=L"$P=0.0,\, \lambda_0=0.001\,\mathrm{eV}$")
ax2.plot(dat3[2:end,1]*fac_t,(dat3[2:end,2]-dat3[2:end,3])*fac_c,color="green",linewidth=lw,label=L"$P=0.5,\, \lambda_0=0.0\,\mathrm{eV}$")
ax2.plot(dat4[2:end,1]*fac_t,(dat4[2:end,2]-dat4[2:end,3])*fac_c,color="orange",linestyle="dashed",linewidth=lw,label=L"$P=0.0,\, \lambda_0=0.0\,\mathrm{eV}$")
ax2.set_ylim(-0.3,0.1)
ax2.set_xlim(1,1.8)

plt.savefig("MF_noPh_2x4_current_dynamics_V0.pdf")

