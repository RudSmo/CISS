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

w0 = 1.0
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

curMFA = readdlm("Asy_g/MF/currents_avg_n30.dat")
popMFA = readdlm("Asy_g/MF/electrons_avg_n30.dat")
curMFS = readdlm("Sym_g/MF/currents_avg_n30.dat")
popMFS = readdlm("Sym_g/MF/electrons_avg_n30.dat")


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

ax1.plot(curMFA[2:end,1]*fac_t, (curMFA[2:end,2]-curMFA[2:end,3])*fac_c*2, color="k",label = "MF-HEOM, anti-sym.", linewidth=lw+2)
ax1.plot(curMFS[2:end,1]*fac_t, (curMFS[2:end,2]-curMFS[2:end,3])*fac_c*2, color="k",linestyle="dashed",label = "MF-HEOM, sym.", linewidth=lw+2)

ax11.plot(curMFA[2:end,1]*fac_t, (curMFA[2:end,2]-curMFA[2:end,3])*fac_c*2, color="k", linewidth=lw+2)
ax11.plot(curMFS[2:end,1]*fac_t, (curMFS[2:end,2]-curMFS[2:end,3])*fac_c*2, color="k", linestyle="dashed", linewidth=lw+2)

ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

ax2.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,2], color="k", linewidth=lw+2)
ax2.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,2], color="k", linestyle="dashed", linewidth=lw+2)
ax22.plot(popMFA[2:end,1]*fac_t, popMFA[2:end,2], color="k", linewidth=lw+2)
ax22.plot(popMFS[2:end,1]*fac_t, popMFS[2:end,2], color="k", linestyle="dashed", linewidth=lw+2)

ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

colors = [("orange","solid"), ("orange","dashed"), ("green","solid"), ("green","dashed")]

for (cur_fil,pop_fil,label,col) in zip(TT_Cur_Files,TT_Pop_Files,TT_labels,colors)
	cur_dat = readdlm(cur_fil)
	pop_dat = readdlm(pop_fil)
	ax1.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3]+cur_dat[2:end,4]-cur_dat[2:end,5])*fac_c,color=col[1],linestyle=col[2],label=label,linewidth=lw)
	ax11.plot(cur_dat[2:end,1]*fac_t, (cur_dat[2:end,2]-cur_dat[2:end,3]+cur_dat[2:end,4]-cur_dat[2:end,5])*fac_c,color=col[1],linestyle=col[2],label=label,linewidth=lw)
	ax2.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col[1],linestyle=col[2],linewidth=lw)
	ax22.plot(pop_dat[2:end,1]*fac_t, pop_dat[2:end,2],color=col[1],linestyle=col[2],linewidth=lw)
end

ax1.legend(loc="best")
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(0,400)
ax1.set_ylim(-10,40)
ax11.set_xlim(100,400)
ax11.set_ylim(27,28.3)
#ax11.set_yticks([13.8,13.9,14.0])

ax2.set_ylabel(L"$n$")
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
#ax2.set_ylim(-1,1)
ax2.set_ylim(-0.1,1.1)
ax2.set_xlim(0,400)
ax22.set_xlim(100,400)
ax22.set_ylim(0.74,0.76)
#ax22.set_yticks([0.748,0.749,0.75])
#ax22.set_yticklabels(["0.748","0.749","0.750"])

plt.savefig("Comparison_TTHEOM_MFHEOM_Holstein_1x1_Spin_w0_$(w0)_V$(V).pdf")

