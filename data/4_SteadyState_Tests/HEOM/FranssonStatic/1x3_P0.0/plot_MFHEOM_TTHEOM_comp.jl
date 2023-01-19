using LinearAlgebra
using PyPlot
using DelimitedFiles
using HDF5

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


fig = plt.figure(figsize = (10,10))
px = 0.102
py = 0.79
sx = 0.88
dy = 0.07
dy1 = dy-0.05
sy = 0.19

fac_c = 1e3 * 2.8494561
fac_t = 1/(4.135668)
fac_s = 1e3

ax1 = plt.axes([px,py,sx,sy+0.007])
ax11 = plt.axes([0.25,0.815,0.4,0.08])
ax1.set_xticklabels([])
ax2 = plt.axes([px,py-sy-dy1,sx,sy+0.007])
ax3 = plt.axes([px,py-2*sy-dy-dy1,sx,sy])
ax4 = plt.axes([px,py-3*sy-2*dy-dy1,sx,sy])


MF_dat_cur = readdlm("MFHEOM/V0.48/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.0_X0_-10.0_np40_dt0.0125_1_currents.dat")
ax1.plot(fac_t * MF_dat_cur[2:end,1],fac_c * (MF_dat_cur[2:end,2]-MF_dat_cur[2:end,3]), color="k", linewidth=lw,label="exact")
ax11.plot(fac_t * MF_dat_cur[2:end,1],fac_c * (MF_dat_cur[2:end,2]-MF_dat_cur[2:end,3]), color="k", linewidth=lw,label="exact")

MF_dat_pop = readdlm("MFHEOM/V0.48/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.0_X0_-10.0_np40_dt0.0125_1_electrons.dat")
ax2.plot(fac_t * MF_dat_pop[2:end,1],fac_s * (MF_dat_pop[2:end,2]-MF_dat_pop[2:end,3]), color="k", linewidth=lw)

ranks = [100,120,140,160,200]
colors = ["red","blue","green","orange","pink"]
for (r,col) in zip(ranks,colors)
	TT_dat_cur = readdlm("np10_r$(r)/V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_$(r).txt")	
	ax1.plot(fac_t * TT_dat_cur[2:end,1],fac_c * (TT_dat_cur[2:end,2]-TT_dat_cur[2:end,3]), color = col, linewidth = lw-1,label=L"$r=\,$"*string(r))
	ax11.plot(fac_t * TT_dat_cur[2:end,1],fac_c * (TT_dat_cur[2:end,2]-TT_dat_cur[2:end,3]), color=col, linewidth=lw-1)

	TT_dat_pop = readdlm("np10_r$(r)/V0.48/n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_$(r).txt") 
	ax2.plot(fac_t * TT_dat_pop[2:end,1],fac_s * (TT_dat_pop[2:end,2]-TT_dat_pop[2:end,3]), color = col, linewidth = lw-1)
end
ax1.legend(loc="lower right",ncol=2)
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
an_x_1 = 0.02
an_y_1 = 0.1
ax1.annotate("(a)",(an_x_1,an_y_1),xycoords="axes fraction",fontsize=fs)
ax11.set_xlim(200,500)
ax11.set_ylim(13.9,14.1)

ax2.set_ylabel(L"$s_1^z\,[\times 10^{-3}]$",labelpad=3)
ax2.set_ylim(-0.5,0.5)
ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax2.annotate("(b)",(an_x_1,an_y_1),xycoords="axes fraction",fontsize=fs)

times = [200,401,1002,1604,4002]
t = 2002
ranks = [100,120,140,160]#,200]
Ir = []
ranks_c = []
for r in ranks
	TT_dat_cur = readdlm("np10_r$(r)/V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_$(r).txt")
	if t < length(TT_dat_cur[2:end,1])
		push!(ranks_c,r)
		push!(Ir,fac_c * (TT_dat_cur[t,2]-TT_dat_cur[t,3]))
	end
end

t_mf = (t-1)*2

ax3.plot(ranks_c,Ir,color="blue",marker="x",markersize=12,label="TT-HEOM")
ax3.plot(ranks,fac_c*(MF_dat_cur[t_mf,2]-MF_dat_cur[t_mf,3])*ones(length(ranks)),color="k",linestyle="dashed",label="exact",linewidth=lw)
ax3.set_xlim(95,165)
ax3.set_xlabel(L"rank $r$")
ax3.legend(loc="upper right")
ytx = [14.70,14.72,14.74,14.76,14.78]
ax3.set_yticks(ytx)
ax3.set_yticklabels(["14.70","14.72","14.74","14.76","14.78"])
ax3.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax3.annotate("(c)",(an_x_1,an_y_1),xycoords="axes fraction",fontsize=fs)

an_x_2 = 0.85
an_y_2 = 0.5
ax3.annotate(L"$t=500\,\mathrm{fs}$",(an_x_2,an_y_2),xycoords="axes fraction",fontsize=fs)
ax3.annotate(L"$N_p=10$",(an_x_2-0.019,an_y_2-0.12),xycoords="axes fraction",fontsize=fs)

np = [10,15,20]
In = []
np_c = []
for n in np
	TT_dat_cur = readdlm("np$(n)_r120/V0.48/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")
	if t < length(TT_dat_cur[2:end,1])
		push!(np_c,n)
		push!(In,fac_c * (TT_dat_cur[t,2]-TT_dat_cur[t,3]))
	end
end
ytx = [14.69,14.70,14.71,14.72]
an_y_2 = 0.85
ax4.plot(np_c,In,color="blue",marker="x",markersize=12)
ax4.plot(np,fac_c*(MF_dat_cur[t_mf,2]-MF_dat_cur[t_mf,3])*ones(length(np)),color="k",linestyle="dashed",linewidth=lw)
ax4.set_xlim(9.5,20.5)
ax4.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax4.set_xlabel(L"Padé poles $N_p$")
ax4.set_yticks(ytx)
ax4.set_yticklabels(["14.69","14.70","14.71","14.72"])
ax4.annotate(L"$t=500\,\mathrm{fs}$",(an_x_2,an_y_2),xycoords="axes fraction",fontsize=fs)
ax4.annotate(L"$r=120$",(an_x_2,an_y_2-0.12),xycoords="axes fraction",fontsize=fs)
ax4.annotate("(d)",(an_x_1,an_y_1),xycoords="axes fraction",fontsize=fs)


ax1.set_xlim(0,500)
ax2.set_xlim(0,500)
		
plt.savefig("Comparison_TTHEOM_MFHEOM_FranssonStatic.pdf")


