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
w0 = 0.05
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

curMFA = readdlm("Asy_g/MF/currents_avg_n30.dat")
curMFS = readdlm("Sym_g/MF/currents_avg_n30.dat")
popMFA = readdlm("Asy_g/MF/electrons_avg_n30.dat")
popMFS = readdlm("Sym_g/MF/electrons_avg_n30.dat")
curNo = readdlm("noPh/L1_gPhU0.0_gPhD0.0_V0.48_w00.002_X0_0.0_dt0.05_1_currents.dat")
popNo = readdlm("noPh/L1_gPhU0.0_gPhD0.0_V0.48_w00.002_X0_0.0_dt0.05_1_electrons.dat")

TT_Cur_Files = ["Asy_g/TT/np10_r100_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
		"Sym_g/TT/np10_r100_nV10_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",	
		"Asy_g/TT/np10_r150_nV15_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
                "Sym_g/TT/np10_r150_nV15_V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt"]
TT_Pop_Files = ["Asy_g/TT/np10_r100_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "Sym_g/TT/np10_r100_nV10_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "Asy_g/TT/np10_r150_nV15_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt",
                "Sym_g/TT/np10_r150_nV15_V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_150.txt"]

TT_labels = [L"$r=100, \ n_\mathrm{vib} = 10$, anti-sym.",
             L"$r=100, \ n_\mathrm{vib} = 10$, sym.",
	     L"$r=150, \ n_\mathrm{vib} = 15$, anti-sym.",
	     L"$r=150, \ n_\mathrm{vib} = 15$, sym."]

ax1.plot(curMFA[2:end,1]*fac_t, (curMFA[2:end,2]-curMFA[2:end,3])*fac_c*2, color="k",label = "MF-HEOM anti-sym.", linewidth=lw+2)
ax1.plot(curMFS[2:end,1]*fac_t, (curMFS[2:end,2]-curMFS[2:end,3])*fac_c*2, color="k",linestyle="dashed",label = "MF-HEOM sym.", linewidth=lw+2)
ax1.plot(curNo[2:end,1]*fac_t, (curNo[2:end,2]-curNo[2:end,3])*fac_c*2, color="orange",label = L"exact $g=0\,\mathrm{meV}$", linewidth=lw,linestyle="dashed")

ax11.plot(curMFA[2:end,1]*fac_t, (curMFA[2:end,2]-curMFA[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+2)
ax11.plot(curNo[2:end,1]*fac_t, (curNo[2:end,2]-curNo[2:end,3])*fac_c*2, color="orange",label = L"exact $g=0\,\mathrm{meV}$", linewidth=lw,linestyle="dashed")
ax11.plot(curMFS[2:end,1]*fac_t, (curMFS[2:end,2]-curMFS[2:end,3])*fac_c*2, color="k",linestyle="dashed",linewidth=lw)

ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

ax2.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,2], color="k", linewidth=lw+1)
ax2.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,2], color="k",linestyle="dashed", linewidth=lw)
ax2.plot(popNo[2:end,1]*fac_t, popNo[2:end,2], color="orange", linewidth=lw,linestyle="dashed")

ax22.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,2], color="k", linewidth=lw+1)
ax22.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,2], color="k",linestyle="dashed", linewidth=lw)
ax22.plot(popNo[2:end,1]*fac_t, popNo[2:end,2], color="orange", linewidth=lw,linestyle="dashed")

ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

ax3.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,3], color="k", linewidth=lw+1)
ax3.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,3], color="k",linestyle="dashed", linewidth=lw)
ax3.plot(popNo[2:end,1]*fac_t, popNo[2:end,2], color="orange", linewidth=lw+1,linestyle="dashed")
ax33.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,3], color="k", linewidth=lw+1)
ax33.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,3], color="k",linestyle="dashed", linewidth=lw)
ax33.plot(popNo[2:end,1]*fac_t, popNo[2:end,2], color="orange", linewidth=lw+1,linestyle="dashed")

ax3.annotate("(c)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)


colors = [["red","solid"],["red","dashed"],["blue","solid"],["blue","dashed"]]
for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)
	pop_dat = readdlm(pop_fil)
	ax1.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3]+cur_dat[2:end,4]-cur_dat[2:end,5])*fac_c,color=col[1],linestyle=col[2],label=label,linewidth=lw)
	ax11.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3]+cur_dat[2:end,4]-cur_dat[2:end,5])*fac_c,color=col[1],linestyle=col[2],label=label,linewidth=lw)
	ax2.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col[1],linestyle=col[2],linewidth=lw)
	ax22.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col[1],linestyle=col[2],linewidth=lw)
	ax3.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,3],color=col[1],linestyle=col[2],linewidth=lw)
	ax33.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,3],color=col[1],linestyle=col[2],linewidth=lw)
end

t_end = 300
ax1.legend(loc="lower left")
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_ylim(-17,40)
ax1.set_xlim(0,t_end)
ax11.set_xlim(50,t_end)
ax11.set_ylim(27,28.5)
#ax11.set_yticks([14.02,14.04,14.06,14.08])
#ax11.set_yticklabels(["14.02","14.04","14.06","14.08"])

ax2.set_ylabel(L"$n_{\uparrow}$")
ax2.set_xticklabels([])
ax2.set_xlim(0,t_end)
ax22.set_xlim(100,t_end)
ax22.set_ylim(0.74,0.76)
ax2.set_ylim(-0.05,1.1)
ax22.set_yticks([0.74,0.75,0.76])


ax3.set_ylabel(L"$n_{\downarrow}$")
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax3.set_xlim(0,t_end)
ax33.set_xlim(100,t_end)
ax33.set_ylim(0.74,0.76)
ax3.set_ylim(-0.05,1.1)
ax33.set_yticks([0.74,0.75,0.76])

plt.savefig("Comparison_TTHEOM_MFHEOM_Holstein_1x1_Spin_w0_$(w0)_V$(V).pdf")

