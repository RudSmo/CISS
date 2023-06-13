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
py = 0.12
sx = 0.87
sy = 0.85
dy = 0.03

fig = plt.figure(figsize=(9,5))
ax1 = plt.axes([px,py,sx,sy])
V = 0.48
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

dat100_2 = readdlm("L1x2/V0.48/TTHEOM/np12_r100_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
dat100_3 = readdlm("L1x3/V0.48/TTHEOM/np12_r100_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")

dat120_2 = readdlm("L1x2/V0.48/TTHEOM/np12_r120_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")
dat120_3 = readdlm("L1x3/V0.48/TTHEOM/np12_r120_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")

dat140_2 = readdlm("L1x2/V0.48/TTHEOM/np12_r140_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_140.txt")
dat140_3 = readdlm("L1x3/V0.48/TTHEOM/np12_r140_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_140.txt")

dat160_2 = readdlm("L1x2/V0.48/TTHEOM/np12_r160_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")
dat160_3 = readdlm("L1x3/V0.48/TTHEOM/np12_r160_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")

dat180_2 = readdlm("L1x2/V0.48/TTHEOM/np12_r180_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")
dat180_3 = readdlm("L1x3/V0.48/TTHEOM/np12_r180_nV30_V0.48_P0.5/norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")


ax1.plot(dat100_2[2:end,1]*fac_t_v2,dat100_2[2:end,2],color="black",label=L"$r=100$",linewidth=lw+2)
ax1.plot(dat100_3[2:end,1]*fac_t_v2,dat100_3[2:end,2],color="black",linestyle="dashed",linewidth=lw+2)

ax1.plot(dat120_2[2:end,1]*fac_t_v2,dat120_2[2:end,2],color="yellow",label=L"$r=120$",linewidth=lw+1)
ax1.plot(dat120_3[2:end,1]*fac_t_v2,dat120_3[2:end,2],color="yellow",linestyle="dashed",linewidth=lw+1)

ax1.plot(dat140_2[2:end,1]*fac_t_v2,dat140_2[2:end,2],color="green",label=L"$r=140$",linewidth=lw+1)
ax1.plot(dat140_3[2:end,1]*fac_t_v2,dat140_3[2:end,2],color="green",linestyle="dashed",linewidth=lw+1)

ax1.plot(dat160_2[2:end,1]*fac_t_v2,dat160_2[2:end,2],color="turquoise",label=L"$r=160$",linewidth=lw+1)
ax1.plot(dat160_3[2:end,1]*fac_t_v2,dat160_3[2:end,2],color="turquoise",linestyle="dashed",linewidth=lw+1)

ax1.plot(dat180_2[2:end,1]*fac_t_v2,dat180_2[2:end,2],color="orange",label=L"$r=180$",linewidth=lw+1)
ax1.plot(dat180_3[2:end,1]*fac_t_v2,dat180_3[2:end,2],color="orange",linestyle="dashed",linewidth=lw+1)

ax1.legend(loc="lower right",ncol=2)
ax1.set_ylabel(L"$\langle\!\langle 1\vert \rho\rangle\!\rangle$")
ax1.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax1.set_ylim(0.98,1.005)
ax1.set_yticks([0.98,0.99,1.00,1.01])
ax1.set_yticklabels([0.98,0.99,"1.00",1.01])
ax1.set_xlim(0,600)

plt.savefig("Norm_TTHEOM_rankVar_LVar.pdf")

