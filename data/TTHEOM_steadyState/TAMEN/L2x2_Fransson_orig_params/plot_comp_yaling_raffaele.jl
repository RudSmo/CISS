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
fac_H = 41.341 
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.09
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

datY = readdlm("current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
datR00 = readdlm("../1x4_sym_np10_r120_nV40_V0.48_P0.5_corrected/Borrelli/curr00.txt")
datR01 = readdlm("../1x4_sym_np10_r120_nV40_V0.48_P0.5_corrected/Borrelli/curr01.txt")
datR10 = readdlm("../1x4_sym_np10_r120_nV40_V0.48_P0.5_corrected/Borrelli/curr10.txt")
datR11 = readdlm("../1x4_sym_np10_r120_nV40_V0.48_P0.5_corrected/Borrelli/curr11.txt")
datE = readdlm("Ehrenfest/currents_avg_n25.dat")

tend = 600
ax.plot(datE[2:end,1]*fac_t, fac_c*2*(datE[2:end,2]-datE[2:end,3]),color="blue",label="Ehrenfest+HEOM",linewidth=lw)
ax.plot(tend-2,fac_c*2*(datE[end,2]-datE[end,3]),color="blue",marker="o",markersize=10)
ax.plot(datY[2:end,1]*fac_t,fac_c*(datY[2:end,2]-datY[2:end,3]+datY[2:end,4]-datY[2:end,5]),color="k",label="Yalings implementation",linewidth=lw)
ax.plot(datR00[2:end,1]*fac_t*fac_H,fac_c*(datR00[2:end,2] + datR10[2:end,2] - datR01[2:end,2] - datR11[2:end,2]),color="red",label="Raffaeles implementation",linewidth=lw)
#ax.legend()
ax.set_xlim(0,tend)
ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")

plt.savefig("Comparison_L1x4_Fransson_Divergence.pdf")
