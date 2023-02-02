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

w0 = 0.05
V = 0.48
num_sites = 4

px = 0.09
py = 0.56
sx = 0.89
sy = 0.42
dy = 0.05

fig = plt.figure(figsize=(10,7))
ax1 = plt.axes([px,py,sx,sy])
#ax11 = plt.axes([px+0.4,py+0.06,0.46,0.15])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
#ax22 = plt.axes([px+0.4,py+0.06-sy-dy,0.46,0.15])


t_end = 150

TT_Cur_Files = ["1x$(num_sites)_sym_np10_r100_nV40_V0.48_P0.5_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
		"1x$(num_sites)_sym_np10_r120_nV40_V0.48_P0.5_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt"]
TT_Pop_Files = ["1x$(num_sites)_sym_np10_r100_nV40_V0.48_P0.5_corrected/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "1x$(num_sites)_sym_np10_r120_nV40_V0.48_P0.5_corrected/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt"]

cur = readdlm("noPh/L4_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat")
pop = readdlm("noPh/L4_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.05_1_electrons.dat")


TT_labels = [L"$r=100$, $n_\mathrm{vib} = 30$",
             L"$r=120$, $n_\mathrm{vib} = 30$"]

colors = ["green", "blue", "red", "turquoise"]
for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)[2:end,:]
	ti_cur = cur_dat[:,1]
	cur_dat = cur_dat[sortperm(ti_cur),:]
	
	pop_dat = readdlm(pop_fil)[2:end,:]
	ti_pop = pop_dat[:,1]	
	pop_dat = pop_dat[sortperm(ti_pop),:]
	ax1.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3]+cur_dat[:,4]-cur_dat[:,5])*fac_c,color=col,label=label,linewidth=lw)
#	ax11.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3]+cur_dat[:,4]-cur_dat[:,5])*fac_c,color=col,label=label,linewidth=lw)
	ax2.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
#	ax22.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
end


ax1.plot(cur[2:end,1]*fac_t, (cur[2:end,2]-cur[2:end,3])*fac_c*2,color="orange",label=L"exact $t_1 = \lambda_1 = 0\, \mathrm{meV}$",linestyle="dashed",linewidth=lw)
ax2.plot(pop[2:end,1]*fac_t, pop[2:end,2], color="orange", linestyle="dashed",linewidth=lw)
ax1.legend(loc="lower left",numpoints=1)
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(0,t_end)
ax1.set_ylim(0,2.7)
#ax11.set_xlim(0,100)
#ax11.set_ylim(1.19,1.21)
#ax11.set_ylim(1.8,2.4)
#ax11.set_yticks([1.16,1.2,1.24])
#ax11.set_yticklabels(["1.16","1.20","1.24"])

ax2.set_ylabel(L"$n_{1\uparrow}$")
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
#ax2.set_ylim(-1,1)
ax2.set_ylim(-0.1,0.9)
ax2.set_xlim(0,t_end)
#ax22.set_xlim(0,100)
#ax22.set_ylim(0.0,0.5)
#ax22.set_yticks([0.748,0.749,0.75])
#ax22.set_yticklabels(["0.748","0.749","0.750"])

plt.savefig("Comparison_TTHEOM_MFHEOM_Fransson_1x$(num_sites)_w0_$(w0)_V$(V).pdf")

