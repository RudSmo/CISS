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


fig = plt.figure(figsize = (10,12))
px = 0.09
py = 0.81
sx = 0.89
dy = 0.011
sy = 0.178
dx = 0.015
sxc = 0.02
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy]) 
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax4 = plt.axes([px,py-3*(sy+dy),sx,sy])
ax5 = plt.axes([px,py-4*(sy+dy),sx,sy])

cp = plt.cm.jet

folds = ["V$(i)" for i in 0:20]
Vs = Array(0:20)*0.1
lw = 3
fac_e = 1
fac_c = 4
fac_t = 1
a = 20000

sz_1 = zeros(length(folds))
electrons_dat = readdlm("../NoPh_1x1/electrons_steadyState.dat")
Vs_s = electrons_dat[:,1]
sz_1s =electrons_dat[:,2]
sz_2 = zeros(length(folds))
sz_2s = electrons_dat[:,3]
cur = zeros(length(folds))
curs = readdlm("../NoPh_1x1/currents_steadyState.dat")[:,2]*fac_c
Xmin = zeros(length(folds))
Xmax = zeros(length(folds))
Pmin = zeros(length(folds))
Pmax = zeros(length(folds))

for (ind,f) in enumerate(folds)
	dat_pop = readdlm("$(f)/electrons_avg_n20.dat")
	sz_1[ind] = dat_pop[end,2]
	sz_2[ind] = dat_pop[end,3]
	
	dat_cur = readdlm("$(f)/currents_avg_n20.dat")
        cur[ind] = fac_c*(dat_cur[end,2]-dat_cur[end,3])

	dat_ph = readdlm("$(f)/phonons_avg_n20.dat")
	Xmin[ind] = minimum(dat_ph[end-a:end,2])
	Xmax[ind] = maximum(dat_ph[end-a:end,2])

	Pmin[ind] = minimum(dat_ph[end-a:end,3])
	Pmax[ind] = maximum(dat_ph[end-a:end,3])
end


ax1.plot(Vs,sz_1,linewidth=lw,color="k")
ax1.plot(Vs_s,sz_1s,linewidth=lw,color="orange")
ax1.set_ylabel(L"$n_\uparrow$")
ax1.set_ylim(0.18,0.5)
ax1.set_xticklabels([])
ax1.annotate("(a)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)


ax2.plot(Vs,sz_2,linewidth=lw,color="k",label="Ehrenfest")
ax2.set_ylim(0.48,0.8)
ax2.plot(Vs_s,sz_2s,linewidth=lw,color="orange",label=L"Steady state $g=0$")
ax2.set_xticklabels([])
ax2.set_ylabel(L"$n_\downarrow$")
ax2.legend(loc="upper right")
ax2.annotate("(b)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)


ax3.plot(Vs,cur,linewidth=lw,color="k")
ax3.set_ylim(-0.05,1.0)
ax3.plot(Vs_s,curs,linewidth=lw,color="orange")
ax3.set_ylabel(L"$I_L - I_R\, [\Gamma_0]$")
ax3.set_xticklabels([])
ax3.annotate("(c)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)

ax4.plot(Vs,Xmin,linewidth=lw,color="blue",label=L"min $(t\in [15000,20000]t_0^{-1})$")
ax4.plot(Vs,Xmax,linewidth=lw,color="red",label=L"max $(t\in [15000,20000]t_0^{-1})$")
ax4.legend(loc="upper right")

ax4.set_xticklabels([])
ax4.annotate("(d)",(0.015,0.85),xycoords="axes fraction",fontsize=fs)
ax4.set_ylabel(L"$X$")

ax5.plot(Vs,Pmin,linewidth=lw,color="blue")
ax5.plot(Vs,Pmax,linewidth=lw,color="red")
ax5.set_xlabel(L"$V\,[t_0]$")
ax5.annotate("(e)",(0.015,0.85),xycoords="axes fraction",fontsize=fs)
ax5.set_ylabel(L"$P$")


plt.savefig("Pops_currents_phonons_steadyState_1x1_VVar.pdf")





