using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 2.8494561
fac_t = 1/(4.135668)
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

w0 = 0.05
V = 0.48

px = 0.09
py = 0.56
sx = 0.89
sy = 0.42
dy = 0.05

fig = plt.figure(figsize=(10,7))
ax1 = plt.axes([px,py,sx,sy])
ax11 = plt.axes([px+0.4,py+0.06,0.46,0.15])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax22 = plt.axes([px+0.4,py+0.06-sy-dy,0.46,0.15])

curMF = readdlm("MF/mid/currents_avg_n30.dat")
popMF = readdlm("MF/mid/electrons_avg_n30.dat")
curNo = readdlm("noPh/L1_gPh0.0_tPh0.0_V0.48_w00.04_X0_0.0_dt0.025_1_currents.dat") 
popNo = readdlm("noPh/L1_gPh0.0_tPh0.0_V0.48_w00.04_X0_0.0_dt0.025_1_electrons.dat")

TT_Cur_Files = ["TT/np10_r100_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV12_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV15_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV17_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt"]

TT_Pop_Files = ["TT/np10_r100_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV12_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV15_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r100_nV17_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt"]
TT_labels = [L"$n_\mathrm{vib} = 10$",
             L"$n_\mathrm{vib} = 12$",
             L"$n_\mathrm{vib} = 15$",
             L"$n_\mathrm{vib} = 17$"]

ax1.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+2)
ax11.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+2)
ax1.plot(curMF[2:end,1]*fac_t, ones(length(curMF[2:end,1]))*(curNo[end,2]-curNo[end,3])*fac_c*2, color="orange",label = L"exact $g=0\,\mathrm{meV}$", linewidth=lw,linestyle="dashed")
ax11.plot(curMF[2:end,1]*fac_t, ones(length(curMF[2:end,1]))*(curNo[end,2]-curNo[end,3])*fac_c*2, color="orange",linewidth=lw,linestyle="dashed")
ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

ax2.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+2)
ax22.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+2)
ax2.plot(popMF[2:end,1]*fac_t, popNo[end,2]*ones(length(popMF[2:end,1])), color="orange", linewidth=lw,linestyle="dashed")
ax22.plot(popMF[2:end,1]*fac_t, popNo[end,2]*ones(length(popMF[2:end,1])), color="orange", linewidth=lw,linestyle="dashed")
ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

colors = ["green", "blue", "red", "turquoise"]
for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)[2:end,:]
	ti_cur = cur_dat[:,1]
	cur_dat = cur_dat[sortperm(ti_cur),:]
	
	pop_dat = readdlm(pop_fil)[2:end,:]
	ti_pop = pop_dat[:,1]	
	pop_dat = pop_dat[sortperm(ti_pop),:]
	ax1.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3])*fac_c,color=col,label=label,linewidth=lw)
	ax11.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3])*fac_c,color=col,label=label,linewidth=lw)
	ax2.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
	ax22.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
end
t_end = 1000
ax1.legend(loc="lower left")
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(0,t_end)
ax1.set_ylim(-10,23)
ax11.set_xlim(100,t_end)
ax11.set_ylim(8,14.5)
#ax11.set_yticks([13.8,13.9,14.0])

ax2.set_ylabel(L"$n$")
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
#ax2.set_ylim(-1,1)
ax2.set_ylim(-0.1,1.1)
ax2.set_xlim(0,t_end)
ax22.set_xlim(100,t_end)
ax22.set_ylim(0.72,0.86)
#ax22.set_yticks([0.748,0.749,0.75])
#ax22.set_yticklabels(["0.748","0.749","0.750"])

plt.savefig("Comparison_TTHEOM_MFHEOM_Holstein_1x1_noSpin_w0_$(w0)_V$(V).pdf")

