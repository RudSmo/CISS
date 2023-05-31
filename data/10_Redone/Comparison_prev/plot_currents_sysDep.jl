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

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.69
sx = 0.87
sy = 0.27
dy = 0.03

fig = plt.figure(figsize=(9,9))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
V = 55
col = ["black","blue","red"]

ws = 0.002
T0 = 0.04/ws

fac_c_v2 = 1e3 * 0.2433946
fac_t_v2 = 1/(4.135668*1000)
fac_c_v3 = 0.002 * 1e3 * 0.2433946
fac_t_v3 = 1/(ws * 4.135668*1000)

Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\mu\mathrm{A}]$",L"$s_1^z$"]

dat13_v2 = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L1x3/V55/currents_avg_n51.dat")
dat13_v3 = readdlm("L1x3/V55/MFHEOM_v3/currents_avg_n41.dat")
dat13_B = readdlm("L1x3/V55/MFHEOM_v3_bundles/currents_avg_n6.dat")

dat24_v2 = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L2x4/V55/currents_avg_n1515_comb.dat")
dat24_v3 = readdlm("L2x4/V55/MFHEOM_v3/P/currents_avg_n41.dat")
dat24_B = readdlm("bundles/V55/L2x4/currents_avg_n5_P.dat")

dat34_v2 = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L3x4/V55/currents_avg_n201.dat")
dat34_B = readdlm("bundles/V55/L3x4/currents_avg_n5_P.dat")

ax1.plot(dat13_v2[2:end,1]*fac_t_v2,(dat13_v2[2:end,2]-dat13_v2[2:end,3])*fac_c_v2,linewidth=lw,color=col[1],label = "v1")
ax1.plot(dat13_v3[2:end,1]*fac_t_v3,(dat13_v3[2:end,2]-dat13_v3[2:end,3])*fac_c_v3,linewidth=lw,color=col[2],label = "v2")
ax1.plot(dat13_B[2:end,1]*fac_t_v3,(dat13_B[2:end,2]-dat13_B[2:end,3])*fac_c_v3,linewidth=lw,color=col[3],label = "bundles")
ax1.legend()
ax1.set_ylabel(ylabels[1])
ax1.set_xticklabels([])
ax1.set_xlim(0,65)
ax1.set_ylim(0.2,1.4)
ax1.annotate(L"$L=1\times 3$",(0.8,0.07),xycoords="axes fraction",fontsize=fs)

ax2.plot(dat24_v2[2:end,1]*fac_t_v2,(dat24_v2[2:end,2]-dat24_v2[2:end,3])*fac_c_v2,linewidth=lw,color=col[1])
ax2.plot(dat24_v3[2:end,1]*fac_t_v3,(dat24_v3[2:end,2]-dat24_v3[2:end,3])*fac_c_v3,linewidth=lw,color=col[2])
ax2.plot(dat24_B[2:end,1]*fac_t_v3,(dat24_B[2:end,2]-dat24_B[2:end,3])*fac_c_v3,linewidth=lw,color=col[3])
ax2.annotate(L"$L=2\times 4$",(0.8,0.07),xycoords="axes fraction",fontsize=fs)
ax2.set_ylabel(ylabels[1])
ax2.set_xlim(0,65)
ax2.set_ylim(0.2,1.4)
ax2.set_xticklabels([])

ax3.plot(dat34_v2[2:end,1]*fac_t_v2,(dat34_v2[2:end,2]-dat34_v2[2:end,3])*fac_c_v2,linewidth=lw,color=col[1])
ax3.plot(dat34_B[2:end,1]*fac_t_v3,(dat34_B[2:end,2]-dat34_B[2:end,3])*fac_c_v3,linewidth=lw,color=col[3])
ax3.set_ylabel(ylabels[1])
ax3.annotate(L"$L=3\times 4$",(0.8,0.07),xycoords="axes fraction",fontsize=fs)
ax3.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax3.set_ylim(0.2,1.4)
ax3.set_xlim(0,65)


plt.savefig("Comparison_MFHEOM_v3_bundles_v2_currents_sysDep.pdf")

