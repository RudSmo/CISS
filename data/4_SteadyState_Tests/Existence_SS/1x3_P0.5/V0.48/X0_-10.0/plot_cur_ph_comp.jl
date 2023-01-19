using LinearAlgebra
using DelimitedFiles
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
lw = 3
al = 0.4
fac_t = 1/(4.135668*1e3)
fac_c = 1e3 * 2.8494561
fac_ph = 0.1973


dat_c = readdlm("L3_tPh0.004_SOC0.001_SOCPh0.0001_V0.48_w00.002_P0.5_X0_-10.0_dt0.025_1_currents.dat")
dat_p = readdlm("L3_tPh0.004_SOC0.001_SOCPh0.0001_V0.48_w00.002_P0.5_X0_-10.0_dt0.025_1_phonons.dat")


fig = plt.figure(figsize=(10,6))
px = 0.1
py = 0.55
sx = 0.89
sy = 0.4
dy = 0.03
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])


ax1.plot(dat_c[2:end,1]*fac_t , (dat_c[2:end,2]-dat_c[2:end,3])*fac_c, linewidth = lw, color = "k")
ax1.set_ylabel(L"$I_L-I_R\,[\mu \mathrm{A}]$")
ax1.set_xticklabels([])
ax1.set_xlim(0,6)
ax1.annotate("(a)",(0.01,0.87),xycoords="axes fraction",fontsize=fs)

ax2.plot(dat_p[2:end,1]*fac_t , dat_p[2:end,2]*fac_ph, linewidth = lw, label = "Bond between sites 1,2")
ax2.plot(dat_p[2:end,1]*fac_t , dat_p[2:end,4]*fac_ph, linewidth = lw, label = "Bond between sites 2,3")
ax2.set_ylabel(L"$X\,[\mu\mathrm{m}]$")
ax2.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax2.annotate("(b)",(0.01,0.87),xycoords="axes fraction",fontsize=fs)

ax2.set_xlim(0,6)

ax2.legend(loc="lower left")

plt.savefig("Fransson_L1x3_V0.48_cur_ph_comp.pdf")
