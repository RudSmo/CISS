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


px = 0.09
py = 0.56
sx = 0.89
sy = 0.42
dy = 0.05

fig = plt.figure(figsize=(10,7))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

curMF = readdlm("MFHEOM/currents_avg_n60.dat")
popMF = readdlm("MFHEOM/electrons_avg_n60.dat")


TT_Cur_Files = ["TTHEOM/np10_r130_nV5_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
		"TTHEOM/np10_r130_nV30_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
		"TTHEOM/np10_r150_nV5_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
		"TTHEOM/np10_r150_nV30_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
		"TTHEOM/np10_r200_nV15_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_200.txt",
 		"TTHEOM/np10_r200_nV30_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_200.txt"]

TT_Pop_Files = ["TTHEOM/np10_r130_nV5_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
                "TTHEOM/np10_r130_nV30_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_130.txt",
                "TTHEOM/np10_r150_nV5_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
                "TTHEOM/np10_r150_nV30_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
                "TTHEOM/np10_r200_nV15_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_200.txt",
                "TTHEOM/np10_r200_nV30_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_200.txt"]

TT_labels = [L"$r=130, \ n_\mathrm{vib} = 5$",
             L"$r=130, \ n_\mathrm{vib} = 30$",
	     L"$r=150, \ n_\mathrm{vib} = 5$",
	     L"$r=150, \ n_\mathrm{vib} = 30$",
	     L"$r=200, \ n_\mathrm{vib} = 15$",
	     L"$r=200, \ n_\mathrm{vib} = 30$"]

ax1.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c, color="k",label = "MF-HEOM", linewidth=lw)
ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

#ax2.plot(popMF[2:end,1]*fac_t, (popMF[2:end,2]-popMF[2:end,3])*0.5, color="k", linewidth=lw)
ax2.plot(popMF[2:end,1]*fac_t, popMF[2:end,3], color="k", linewidth=lw)
ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

colors = ["green", "blue", "red", "turquoise", "orange", "purple"]
for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)
	pop_dat = readdlm(pop_fil)
	ax1.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]+cur_dat[2:end,4]-cur_dat[2:end,3]-cur_dat[2:end,5])*fac_c,color=col,label=label,linewidth=lw)
	ax2.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,3],color=col,linewidth=lw)
end

ax1.legend(loc="best")
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(0,50)
ax1.set_ylim(-300,300)

ax2.set_ylabel(L"$n_{1\downarrow}$")
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
#ax2.set_ylim(-1,1)
ax2.set_ylim(0.0,3.5)
ax2.set_xlim(0,50)

plt.savefig("Comparison_TTHEOM_MFHEOM_1x3_w0_0.05_Gamma0.1t0.pdf")

