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
px = 0.08
py = 0.71
sx1 = 0.4
dx1 = 0.09
sx = 0.43
dx = 0.034
dy = 0.05
sy = 0.275
ax11 = plt.axes([px,py,sx1,sy])
ax12 = plt.axes([px+sx1+dx1,py,sx1,sy]) 
ax2 = plt.axes([px,py-sy-dy,2*sx+dx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),2*sx+dx,sy])

folds = ["V0","V4","V12","V20"]#,"V16","V20"]
foldsSS = [0.0,0.4,1.2,2.0]#,1.6,2.0]
fold_str = [L"$0.0$",L"$0.4$",L"$1.2$",L"$2.0$"]#,L"$1.6$",L"$2.0$"]
lw = 2
fac_e = 1
fac_c = 4
fac_t = 1

#cols = plt.cm.Set1(0:1/(length(folds)-1):1)
#cols = [cols[i,:] for i in 1:length(cols[:,1])]
cols = ["red","green","blue","orange"]

ms = 10
for (f,fstr,col) in zip(folds,fold_str,cols)
	
	dat_pop = readdlm("$(f)/electrons_avg_n10.dat")
	ax11.plot(dat_pop[2:end,1]*fac_t,fac_e*(dat_pop[2:end,2]-dat_pop[2:end,3]),linewidth=lw,color=col)
	#ax11.plot(dat_pop[2:end,1]*fac_t,fac_e*dat_pop[2:end,3],linewidth=lw,linestyle="dashed",color=col)
	ax12.plot(dat_pop[2:end,1]*fac_t,fac_e*(dat_pop[2:end,4]-dat_pop[2:end,5]),linewidth=lw,label=L"$V=\,$"*fstr,color=col)
	#ax12.plot(dat_pop[2:end,1]*fac_t,fac_e*dat_pop[2:end,5],linewidth=lw,linestyle="dashed",color=col)

	dat_cur = readdlm("$(f)/currents_avg_n10.dat")
        ax2.plot(dat_cur[2:end,1]*fac_t,(dat_cur[2:end,2]-dat_cur[2:end,3])*fac_c,linewidth=lw,color=col)

	

	dat_ph = readdlm("$(f)/phonons_avg_n10.dat")
	ax3.plot(dat_ph[2:end,2],dat_ph[2:end,3],linewidth=lw,color=col)
	ax3.plot(dat_ph[2,2],dat_ph[2,3],marker="o",markersize=ms,color=col)
end


plt_ti = 1980
for (fSS,col) in zip(foldsSS,cols)
	
	datSS_E = readdlm("NoPh_1x2/V$(fSS)/L2_tPh0.0_SOC0.0_SOCPh0.0_V$(fSS)_w00.05_P0.0_X0_0.0_dt0.002_1_electrons.dat")
	ax11.plot([plt_ti],datSS_E[end,2]-datSS_E[end,3],marker="o",markersize=ms,color=col)
	ax12.plot([plt_ti],datSS_E[end,4]-datSS_E[end,5],marker="o",markersize=ms,color=col)	
	
	datSS_C = readdlm("NoPh_1x2/V$(fSS)/L2_tPh0.0_SOC0.0_SOCPh0.0_V$(fSS)_w00.05_P0.0_X0_0.0_dt0.002_1_currents.dat")
	ax2.plot([plt_ti],(datSS_C[end,2]-datSS_C[end,3])*fac_c,marker="o",markersize=ms,color=col)

end

lw = 2



ax11.set_ylabel(L"$s^z$")
ax11.set_ylim(-0.5,0.05)
xlim_list = [0,500,1000,1500,2000]
xtick_list = [0.0,0.5,1.0,1.5,2.0]
ax11.set_xticks(xlim_list)
ax11.set_xticklabels(xtick_list)
ax11.set_xlabel(L"$t\, [\times 10^3 t_0^{-1}]$",fontsize=13,labelpad=-2)
ax11.annotate("(a)",(0.04,0.9),xycoords="axes fraction",fontsize=fs)
ax11.annotate("site 1",(0.8,0.9),xycoords="axes fraction",fontsize=fs)

#ax12.set_ylabel(L"$n_2$")
ax12.set_ylim(-0.5,0.05)
ax12.set_xticks(xlim_list)
ax12.set_xticklabels(xtick_list)
ax12.set_xlabel(L"$t\, [\times 10^3 t_0^{-1}]$",fontsize=13,labelpad=-2)
ax12.annotate("(b)",(0.04,0.9),xycoords="axes fraction",fontsize=fs)
ax12.annotate("site 2",(0.8,0.9),xycoords="axes fraction",fontsize=fs)
ax12.legend(loc="lower left",numpoints=1,ncol=2)

ax2.set_xticks(xlim_list)
ax2.set_xticklabels(xtick_list)
ax2.set_ylabel(L"$I_L-I_R \, [\Gamma_0]$",labelpad=-2)
ax2.set_xlabel(L"$t\, [\times 10^3 t_0^{-1}]$",fontsize=15,labelpad=-2)
ax2.annotate("(c)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)

ax3.annotate("(d)",(0.015,0.85),xycoords="axes fraction",fontsize=fs)
ax3.set_ylabel(L"$P$")
ax3.set_xlabel(L"$X$")


plt.savefig("Pops_currents_phonons_Peierls_1x2_VVar.pdf")





