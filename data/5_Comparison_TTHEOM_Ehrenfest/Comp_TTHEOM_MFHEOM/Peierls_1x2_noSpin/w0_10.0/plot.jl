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
rc("legend",fontsize=fs-3)
w0 = 10.0
V = 0.48

px = 0.09
py = 0.69
sx = 0.89
sy = 0.2972
dy = 0.012

fig = plt.figure(figsize=(10,9))
ax1 = plt.axes([px,py,sx,sy])
ax11 = plt.axes([px+0.4,py+0.05,0.45,0.12])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax22 = plt.axes([px+0.4,py+0.05-sy-dy,0.45,0.12])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax33 = plt.axes([px+0.4,py+0.05-2*(sy+dy),0.45,0.12])

curMF = readdlm("MF/currents_avg_n30.dat")
popMF = readdlm("MF/electrons_avg_n30.dat")


TT_Cur_Files = ["TT/np10_r100_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
		"TT/np10_r130_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
		"TT/np10_r150_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
		"TT/np10_r150_nV15_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt"]

TT_Pop_Files = ["TT/np10_r100_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TT/np10_r130_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
                "TT/np10_r150_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
                "TT/np10_r150_nV15_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt"]

TT_labels = [L"$r=100, \ n_\mathrm{vib} = 10$",
             L"$r=130, \ n_\mathrm{vib} = 10$",
	     L"$r=150, \ n_\mathrm{vib} = 10$",
	     L"$r=150, \ n_\mathrm{vib} = 15$"]

ax1.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+1)
ax11.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+1)
ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

#ax2.plot(popMF[2:end,1]*fac_t, (popMF[2:end,2]-popMF[2:end,3])*0.5, color="k", linewidth=lw)
ax2.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+1)
ax22.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+1)
ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

ax3.plot(popMF[2:end,1]*fac_t, popMF[2:end,3], color="k", linewidth=lw+1)
ax33.plot(popMF[2:end,1]*fac_t, popMF[2:end,3], color="k", linewidth=lw+1)
ax3.annotate("(c)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)


colors = ["green", "blue", "red", "turquoise"]
for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)
	pop_dat = readdlm(pop_fil)
	ax1.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3])*fac_c,color=col,label=label,linewidth=lw)
	ax11.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3])*fac_c,color=col,label=label,linewidth=lw)
	ax2.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col,linewidth=lw)
	ax22.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col,linewidth=lw)
	ax3.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,3],color=col,linewidth=lw)
	ax33.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,3],color=col,linewidth=lw)
end

ax1.legend(loc="lower left")
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_ylim(-6,25)
ax11.set_xlim(100,250)
ax11.set_ylim(14,14.05)
ax11.set_yticks([14.0,14.02,14.04])
ax11.set_yticklabels(["14.00","14.02","14.04"])

ax2.set_ylabel(L"$n_{1}$")
ax2.set_xticklabels([])
ax22.set_xlim(100,250)
ax22.set_ylim(0.74,0.76)
ax2.set_ylim(-0.05,1.1)
ax22.set_yticks([0.74,0.75,0.76])


ax3.set_ylabel(L"$n_{2}$")
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax33.set_xlim(100,250)
ax33.set_ylim(0.73,0.75)
ax3.set_ylim(-0.05,1.1)
ax33.set_yticks([0.73,0.74,0.75])

plt.savefig("Comparison_TTHEOM_MFHEOM_Peierls_1x2_noSpin_w0_$(w0)_V$(V).pdf")

