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
P = 0.5
px = 0.09
py = 0.69
sx = 0.888
sy = 0.295
dy = 0.02

fig = plt.figure(figsize=(10,10))
ax1 = plt.axes([px,py,sx,sy])
ax11 = plt.axes([px+0.4,py+0.06,0.46,0.15])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax22 = plt.axes([px+0.4,py+0.04-sy-dy,0.46,0.15])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])

t_end = 750

#curMF = readdlm("MFHEOM/P$(P)/currents_avg_n30.dat")
#popMF = readdlm("MFHEOM/P$(P)/electrons_avg_n30.dat")
#curNo = readdlm("noPh/L2_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P$(P)_X0_0.0_dt0.05_1_currents.dat") 
#popNo = readdlm("noPh/L2_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P$(P)_X0_0.0_dt0.05_1_electrons.dat")
#curSS = readdlm("Christoph_HEOM/CU2352.dat")
#popSS = readdlm("Christoph_HEOM/EL2352.dat") 
#ax1.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+2)
#ax1.plot(t_end-10,(curSS[end,2]-curSS[end,3])*fac_c_SS,color="orange",marker="o",label="HEOM steady-state")
#ax11.plot(curMF[2:end,1]*fac_t, (curMF[2:end,2]-curMF[2:end,3])*fac_c*2, color="k",label = "MF-HEOM", linewidth=lw+2)
#ax11.plot(t_end-10,(curSS[end,2]-curSS[end,3])*fac_c_SS,color="orange",marker="o")
#ax1.plot(curMF[2:end,1]*fac_t, ones(length(curMF[2:end,1]))*(curNo[end,2]-curNo[end,3])*fac_c*2, color="orange",label = L"exact $t_1=\lambda_1=0\,\mathrm{meV}$", linewidth=lw,linestyle="dashed")
#ax11.plot(curMF[2:end,1]*fac_t, ones(length(curMF[2:end,1]))*(curNo[end,2]-curNo[end,3])*fac_c*2, color="orange",linewidth=lw,linestyle="dashed")

#ax2.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+2)
#ax22.plot(popMF[2:end,1]*fac_t, popMF[2:end,2], color="k", linewidth=lw+2)
#ax2.plot(popMF[2:end,1]*fac_t, popNo[end,2]*ones(length(popMF[2:end,1])), color="orange", linewidth=lw,linestyle="dashed")
#ax22.plot(popMF[2:end,1]*fac_t, popNo[end,2]*ones(length(popMF[2:end,1])), color="orange", linewidth=lw,linestyle="dashed")
#ax2.plot(t_end-10,popSS[end,2],color="orange",marker="o")
#ax22.plot(t_end-10,popSS[end,2],color="orange",marker="o")


curname = "current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_"
popname = "n_site_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_"
normname = "norm_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_"

conv_params = [("10","100"),("15","100"),("20","100"),("20","110"),("20","120"),("20","130"),("20","140"),("30","100")]
dir_pre = ["nV10","nV15","nV20_r100","nV20_r110","nV20_r120","nV20_r130","nV20_r140","nV30"]

ax1.annotate("(a)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)
ax2.annotate("(b)",(0.02,0.9),xycoords="axes fraction",fontsize=fs)

colors = ["green", "blue", "red","black","orange","purple","magenta","gray"]
for (par,dir,col) in zip(conv_params,dir_pre,colors) 
	cur_dat = readdlm("$(dir)/$(curname)$(par[2]).txt")[2:end,:]
	ti_cur = cur_dat[:,1]
	cur_dat = cur_dat[sortperm(ti_cur),:]
	
	pop_dat = readdlm("$(dir)/$(popname)$(par[2]).txt")[2:end,:]
	norm_dat = readdlm("$(dir)/$(normname)$(par[2]).txt")[2:end,:]

	ti_pop = pop_dat[:,1]	
	pop_dat = pop_dat[sortperm(ti_pop),:]
	label = L"$n_v =\, $"*par[1]*L"$, \ r=\, $"*par[2]

	ax1.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3]+cur_dat[:,4]-cur_dat[:,5])*fac_c,color=col,label=label,linewidth=lw)
	ax11.plot(cur_dat[:,1]*fac_t, (cur_dat[:,2]-cur_dat[:,3]+cur_dat[:,4]-cur_dat[:,5])*fac_c,color=col,label=label,linewidth=lw)
	ax2.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
	ax22.plot(pop_dat[:,1]*fac_t, pop_dat[:,2],color=col,linewidth=lw)
	ax3.plot(norm_dat[:,1]*fac_t,norm_dat[:,2],color=col,linewidth=lw)
end


ax1.legend(loc="lower left",numpoints=1)
ax1.set_xticklabels([])
ax1.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(0,t_end)
ax1.set_ylim(-0.5,2.5)
ax11.set_xlim(100,t_end)
ax11.set_ylim(2.12,2.18)
ax11.set_yticks([2.12,2.14,2.16,2.18])
ax11.set_yticklabels(["2.12","2.14","2.16","2.18"])

ax2.set_ylabel(L"$n_{1\uparrow}$")
#ax2.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax2.set_xticklabels([])
#ax2.set_ylim(-1,1)
ax2.set_ylim(-0.1,1.0)
ax2.set_xlim(0,t_end)
ax22.set_xlim(100,t_end)
ax22.set_ylim(0.8,0.815)
ax22.set_yticks([0.8,0.805,0.81,0.815])
#ax22.set_yticklabels(["0.800","0.805","0.810","0.815"])

ax3.set_ylabel(L"Norm $\vert\!\vert \rho\vert\!\vert$")
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax3.set_ylim(0.996,1.006)
ax3.set_yticks([0.996,0.998,1.0,1.002,1.004,1.006])
ax3.set_yticklabels(["0.996","0.998","1.000","1.002","1.004","1.006"])
ax3.set_xlim(0,t_end)
ax3.annotate("(c)",(0.02,0.87),xycoords="axes fraction",fontsize=fs)

plt.savefig("Convergence_TTHEOM_Fransson_1x2_w0_$(w0)_V$(V)_P$(P)_norm.pdf")

