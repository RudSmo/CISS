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
py = 0.55
sx = 0.87
sy = 0.43
dy = 0.03

fig = plt.figure(figsize=(9,7))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
V = 0.48
axes = [ax1,ax2]
colors = ["black","blue","red"]

ws = 0.002
T0 = 0.04/ws

fac_c_v2 = 1e3 * 0.2433946
fac_t_v2 = 1/(4.135668)
fac_c_v3 = 0.002 * 1e3 * 0.2433946
fac_t_v3 = 1/(ws * 4.135668)

Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\mu\mathrm{A}]$",L"$s_1^z$"]


dat100C = readdlm("TTHEOM/np12_r100_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
dat100E = readdlm("TTHEOM/np12_r100_nV30_V0.48_P0.5/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")

dat120C = readdlm("TTHEOM/np12_r120_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")
dat120E = readdlm("TTHEOM/np12_r120_nV30_V0.48_P0.5/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")

dat160C = readdlm("TTHEOM/np12_r160_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")
dat160E = readdlm("TTHEOM/np12_r160_nV30_V0.48_P0.5/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")

dat180C = readdlm("TTHEOM/np12_r180_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")
dat180E = readdlm("TTHEOM/np12_r180_nV30_V0.48_P0.5/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")

dat2C = readdlm("noPh/L3_tPh_0.0_thPd_0.0_SOC0.5_SOCPh0.0_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_currents.dat")
dat2E = readdlm("noPh/L3_tPh_0.0_thPd_0.0_SOC0.5_SOCPh0.0_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_electrons.dat")

dat3C = readdlm("MFHEOM_v3_PVar_InitGS/currents_avg_n3_P.dat")
dat3E = readdlm("MFHEOM_v3_PVar_InitGS/electrons_avg_n3_P.dat")

ax1.plot(dat3C[2:end,1]*fac_t_v3, (dat3C[2:end,2]-dat3C[2:end,3])*fac_c_v3, color=colors[3], label = "Ehrenfest", linewidth=lw+2)
ax1.plot(dat100C[2:end,1]*fac_t_v2, (dat100C[2:end,2]-dat100C[2:end,3]+dat100C[2:end,4]-dat100C[2:end,5])*fac_c_v2/2, color=colors[1], label = "TTHEOM, "*L"$r=100$", linewidth=lw+1)
ax1.plot(dat120C[2:end,1]*fac_t_v2, (dat120C[2:end,2]-dat120C[2:end,3]+dat120C[2:end,4]-dat120C[2:end,5])*fac_c_v2/2, color="yellow", label = "TTHEOM, "*L"$r=120$", linewidth=lw)
ax1.plot(dat160C[2:end,1]*fac_t_v2, (dat160C[2:end,2]-dat160C[2:end,3]+dat160C[2:end,4]-dat160C[2:end,5])*fac_c_v2/2, color="turquoise", label = "TTHEOM, "*L"$r=160$", linewidth=lw)
ax1.plot(dat180C[2:end,1]*fac_t_v2, (dat180C[2:end,2]-dat180C[2:end,3]+dat180C[2:end,4]-dat180C[2:end,5])*fac_c_v2/2, color="orange", label = "TTHEOM, "*L"$r=180$", linewidth=lw)

ax1.plot(dat2C[2:end,1]*fac_t_v3, (dat2C[2:end,2]-dat2C[2:end,3])*fac_c_v3, color=colors[2], label = L"$t_1=\lambda_1=0$", linewidth=lw-1,linestyle="dashed")

ax1.legend(loc="lower right",ncol=2)
ax1.set_ylabel(ylabels[1])
ax1.set_xticklabels([])
ax1.set_ylim(0.9,1.3)
ax1.set_xlim(0,600)
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)



ax2.plot(dat3E[2:end,1]*fac_t_v3, (dat3E[2:end,2]-dat3E[2:end,3]), color=colors[3], label = "Bundles Ehrenfest "*L"$(n=500000)$", linewidth=lw+2)
ax2.plot(dat100E[2:end,1]*fac_t_v2, (dat100E[2:end,2]-dat100E[2:end,3]), color="black", label = "TTHEOM", linewidth=lw)
ax2.plot(dat120E[2:end,1]*fac_t_v2, (dat120E[2:end,2]-dat120E[2:end,3]), color="yellow", label = "TTHEOM", linewidth=lw)
ax2.plot(dat160E[2:end,1]*fac_t_v2, (dat160E[2:end,2]-dat160E[2:end,3]), color="turquoise", label = "TTHEOM", linewidth=lw)
ax2.plot(dat180E[2:end,1]*fac_t_v2, (dat180E[2:end,2]-dat180E[2:end,3]), color="orange", label = "TTHEOM", linewidth=lw)
ax2.plot(dat2E[2:end,1]*fac_t_v3, (dat2E[2:end,2]-dat2E[2:end,3]), color=colors[2], linewidth=lw-1,linestyle="dashed")
ax2.set_ylabel(ylabels[2],labelpad=3)
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax2.set_xlim(0,600)
ax2.annotate("(b)",(annx,anny),xycoords="axes fraction",fontsize=fs)

plt.savefig("L3_comparison_TTHEOM_MFHEOM_v3_PVar.pdf")

