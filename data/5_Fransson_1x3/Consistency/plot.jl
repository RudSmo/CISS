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

px = 0.11
py = 0.12
sx = 0.86
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

w = 0.01
t1 = 0.02

datS = readdlm("Spin/L3_tPh_-0.02_thPd_-0.02_SOC0.0_SOCPh0.0_V0.4_w00.01_P0.0_noEps_dt0.05_1_currents.dat")
datN = readdlm("noSpin/L3_gPh0.0_tPh-0.02_V0.4_w00.01_E-0.24_dt0.05_1_currents.dat")
ax.plot(datS[2:end,1]*fac_t,(datS[2:end,2]-datS[2:end,3])*fac_c/6,color="k",linewidth=lw,label="Spinfull")
ax.plot(datN[2:end,1]*fac_t,(datN[2:end,2]-datN[2:end,3])*fac_c/3,color="red",linewidth=lw,label="Spinless")

#datS = readdlm("Spin/L3_tPh_-0.02_thPd_-0.02_SOC0.0_SOCPh0.0_V0.4_w00.01_P0.0_noEps_dt0.05_1_phonons.dat")
#datN = readdlm("noSpin/L3_gPh0.0_tPh-0.02_V0.4_w00.01_E-0.24_dt0.05_1_phonons.dat")
#
#
#ax.plot(datS[2:end,1]*fac_t,(datS[2:end,2]-datS[2:end,3])/6,color="k",linewidth=lw,label="Spinfull")
#ax.plot(datN[2:end,1]*fac_t,(datN[2:end,2]-datN[2:end,3])/3,color="red",linewidth=lw,label="Spinless")



ax.legend()
ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
ax.set_ylabel(L"$I/N\,[\mu\mathrm{A}]$")
plt.savefig("Peierls_L3_spin_noSpin_comparison.pdf")

#ax.set_ylabel(L"$X_1/N$")
#plt.savefig("Peierls_L3_spin_noSpin_comparison_phonons.pdf")
