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

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

w = 0.01
t1 = 0.02
datTTe = readdlm("TTHEOM/empty/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
datTTf = readdlm("TTHEOM/full/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
ax.plot(datTTe[2:end,1]*fac_t , (datTTe[2:end,2]-datTTe[2:end,3]+datTTe[2:end,4]-datTTe[2:end,5])*fac_c, color = "k", label="TTHEOM "*L"$n_{i\sigma}(t_0)=0$", linewidth=lw-1)
ax.plot(datTTf[2:end,1]*fac_t , (datTTf[2:end,2]-datTTf[2:end,3]+datTTf[2:end,4]-datTTf[2:end,5])*fac_c, color = "k", label="TTHEOM "*L"$n_{i\sigma}(t_0)=1$", linestyle="dashed", linewidth=lw)

for i in 1:15
	dat = readdlm("MFHEOM/data$(i)/L2_tPhU_-0.04_tPhD_-0.04_SOC0.0_SOCPh0.0_V0.48_dt0.05_1_currents.dat")
	if i == 1
		ax.plot(dat[2:end,1]*fac_t,(dat[2:end,2]-dat[2:end,3])*fac_c*2,linewidth=lw-2,alpha=0.3,color="blue",label="MFHEOM trajectories")
	else
		ax.plot(dat[2:end,1]*fac_t,(dat[2:end,2]-dat[2:end,3])*fac_c*2,linewidth=lw-2,alpha=0.3,color="blue")
	end
end

datMFAvg = readdlm("MFHEOM/currents_avg_n15.dat")
ax.plot(datMFAvg[2:end,1]*fac_t,(datMFAvg[2:end,2]-datMFAvg[2:end,3])*fac_c*2,linewidth=lw,color="red",label="MFHEOM averaged")
ax.legend()
ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
ax.set_ylabel(L"$I \,[\mu\mathrm{A}]$")
ax.set_ylim(0,4.)
ax.set_xscale("log")
plt.savefig("Fransson_L2_P0_non-unique_steady_state_currents_MFHEOM_vs_TTHEOM.pdf")


