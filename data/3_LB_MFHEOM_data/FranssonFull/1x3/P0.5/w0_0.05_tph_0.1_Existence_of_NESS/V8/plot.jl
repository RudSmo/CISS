using LinearAlgebra
using PyPlot
using DelimitedFiles

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
py = 0.72
sx1 = 0.285
dx1 = 0.02
sx = 0.43
dx = 0.038
dy = 0.05
sy = 0.275
ax11 = plt.axes([px,py,sx1,sy])
ax12 = plt.axes([px+sx1+dx1,py,sx1,sy]) 
ax13 = plt.axes([px+2*(sx1+dx1),py,sx1,sy])
ax2 = plt.axes([px,py-sy-dy,2*sx+dx,sy])
ax31 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax32 = plt.axes([px+sx+dx,py-2*(sy+dy),sx,sy])

folds = ["X0_0.0","X0_-0.5","X0_-4.0","X0_-12.0","X0_-25.0"]
fold_str = [L"$0$",L"$-0.5$",L"$-4$",L"$-12$",L"$-25$"]
a = 10000
lw = 2
fac = 1e0
fac_c = 1e0
for (f,fstr) in zip(folds,fold_str)
	dat_pop = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V8.0_$(f)_dt0.002_electrons.dat")
	ax11.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,2]-dat_pop[2:end-a,3])*fac,linewidth=lw)
	ax12.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,4]-dat_pop[2:end-a,5])*fac,linewidth=lw)
	ax13.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,6]-dat_pop[2:end-a,7])*fac,linewidth=lw)
	dat_cur = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V8.0_$(f)_dt0.002_currents.dat")
        ax2.plot(dat_cur[2:end,1],(dat_cur[2:end,2]-dat_cur[2:end,3])*fac_c,linewidth=lw,label=L"$X_0=$"*fstr)

	dat_ph = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V8.0_$(f)_dt0.002_phonons.dat")
	ax31.plot(dat_ph[2:end,2],dat_ph[2:end,3],linewidth=lw)
	ax32.plot(dat_ph[2:end,4],dat_ph[2:end,5],linewidth=lw)
end

ax11.set_ylabel(L"$s^z$")
ax11.set_ylim(-0.02,0.18)
ax11.set_xticks([0,1000,2000,3000,4000])
ax11.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax11.annotate("(a)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax11.annotate("site 1",(0.8,0.82),xycoords="axes fraction",fontsize=fs)

ax12.set_ylim(-0.02,0.18)
ax12.set_yticklabels([])
ax12.set_xticks([0,1000,2000,3000,4000])
ax12.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax12.annotate("(b)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax12.annotate("site 2",(0.8,0.82),xycoords="axes fraction",fontsize=fs)

ax13.set_ylim(-0.02,0.18)
ax13.set_yticklabels([])
ax13.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax13.annotate("(c)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax13.annotate("site 3",(0.8,0.82),xycoords="axes fraction",fontsize=fs)

ax2.set_ylim(-0.01,0.1)
ax2.set_ylabel(L"$I_L-I_R \, [t_0]$")
ax2.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=15,labelpad=-2)
ax2.annotate("(d)",(0.02,0.85),xycoords="axes fraction",fontsize=fs)
ax2.legend(loc="upper right",ncol=2)

ax31.annotate("(e)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax31.set_ylabel(L"$P_1$")
ax31.set_xlabel(L"$X_1$")
ax31.set_ylim(-12,12)
ax31.set_xlim(-50,40)

ax32.annotate("(f)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax32.set_yticklabels([])
ax32.set_ylim(-12,12)
ax32.set_xlim(-50,40)
ax32.set_xlabel(L"$X_2$")
ax32.set_ylabel(L"$P_2$",labelpad=-2)

plt.savefig("Pops_currents_phonons_1x3_V8.pdf")





