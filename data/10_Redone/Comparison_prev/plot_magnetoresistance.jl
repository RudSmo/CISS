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
sy = 0.83
dy = 0.03

fig = plt.figure(figsize=(9,5))
ax1 = plt.axes([px,py,sx,sy])
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

function MR(IP,IM)
	return 100 * (IP .- IM) ./ (IP .+ IM)
end

dat_v2_P = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L2x4/V55/currents_avg_n101.dat")
dat_v2_M = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P-0.5/L2x4/V55/currents_avg_n101.dat")

dat_v3_P = readdlm("L2x4/V55/MFHEOM_v3/P/currents_avg_n41.dat")
dat_v3_M = readdlm("L2x4/V55/MFHEOM_v3/M/currents_avg_n41.dat")

dat_B_P = readdlm("bundles/V55/L3x4/currents_avg_n5_P.dat")
dat_B_M = readdlm("bundles/V55/L3x4/currents_avg_n5_M.dat")


ax1.plot(dat_v2_P[2:end,1]*fac_t_v2,MR(dat_v2_P[2:end,2]-dat_v2_P[2:end,3],dat_v2_M[2:end,2]-dat_v2_M[2:end,3]),linewidth=lw,color=col[1],label = "v1")
ax1.plot(dat_v3_P[2:end,1]*fac_t_v3,MR(dat_v3_P[2:end,2]-dat_v3_P[2:end,3],dat_v3_M[2:end,2]-dat_v3_M[2:end,3]),linewidth=lw,color=col[2],label = "v2")
ax1.plot(dat_B_P[2:end,1]*fac_t_v3,MR(dat_B_P[2:end,2]-dat_B_P[2:end,3],dat_B_M[2:end,2]-dat_B_M[2:end,3]),linewidth=lw,color=col[3],label = "bundles")

ax1.legend()
ax1.set_ylabel(L"$M_R\,[\%]$")
#ax1.set_xlim(0,65)
ax1.set_ylim(-10,12)
ax1.set_xlabel(L"$t\,[\mathrm{ps}]$")


plt.savefig("Comparison_MFHEOM_v3_bundles_v2_magnetoresistance_L2x4_V55.pdf")

