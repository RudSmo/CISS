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
px = 0.065
py = 0.81
sx = 0.8
dy = 0.011
sy = 0.178
dx = 0.015
sxc = 0.02
ax1 = plt.axes([px,py,sx,sy])
ax1c = plt.axes([px+sx+dx,py,sxc,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy]) 
ax2c = plt.axes([px+sx+dx,py-sy-dy,sxc,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax3c = plt.axes([px+sx+dx,py-2*(sy+dy),sxc,sy])
ax4 = plt.axes([px,py-3*(sy+dy),sx,sy])
ax4c = plt.axes([px+sx+dx,py-3*(sy+dy),sxc,sy])
ax5 = plt.axes([px,py-4*(sy+dy),sx,sy])
ax5c = plt.axes([px+sx+dx,py-4*(sy+dy),sxc,sy])

cp = plt.cm.jet

folds = ["V$(i)" for i in 0:2:20]
lw = 2
fac_e = 1
fac_c = 4
fac_t = 1

ms = 10
dat = readdlm("V0/currents_avg_n10.dat")
sz_1 = zeros(length(folds),length(dat[2:2:end,1]))
sz_2 = zeros(length(folds),length(dat[2:2:end,1]))
cur = zeros(length(folds),length(dat[2:2:end,1]))
X = zeros(length(folds),length(dat[2:2:end,1]))
P = zeros(length(folds),length(dat[2:2:end,1]))

for (ind,f) in enumerate(folds)
	dat_pop = readdlm("$(f)/electrons_avg_n10.dat")
	sz_1[end+1-ind,:] = dat_pop[2:2:end,2]-dat_pop[2:2:end,3]
	sz_2[end+1-ind,:] = dat_pop[2:2:end,4]-dat_pop[2:2:end,5]
	
	dat_cur = readdlm("$(f)/currents_avg_n10.dat")
        cur[end+1-ind,:] = fac_c*(dat_cur[2:2:end,2]-dat_cur[2:2:end,3])

	dat_ph = readdlm("$(f)/phonons_avg_n10.dat")
	X[end+1-ind,:] = dat_ph[2:2:end,2]
	P[end+1-ind,:] = dat_ph[2:2:end,3]
end


ax1.imshow(sz_1,aspect="auto",extent=[0,dat[end,1],0,2],cmap="jet")
ax1.set_ylabel(L"$V$")
ax1.set_xticklabels([])
ax1.annotate("(a)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)

cmap = plt.cm.ScalarMappable(cmap=cp,norm=plt.Normalize(vmin = minimum(sz_1), vmax = maximum(sz_1)))
cbar1 = fig.colorbar(cmap, ax=ax1, cax = ax1c, aspect=90,pad=.001)
cbar1.set_label(L"$s_1^z$")
cbar1.ax.tick_params(labelsize=fs)


ax2.imshow(sz_2,aspect="auto",extent=[0,dat[end,1],0,2],cmap="jet")
ax2.set_xticklabels([])
ax2.set_ylabel(L"$V$")
ax2.annotate("(b)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)
cmap = plt.cm.ScalarMappable(cmap=cp,norm=plt.Normalize(vmin = minimum(sz_2), vmax = maximum(sz_2)))
cbar2 = fig.colorbar(cmap, ax=ax2, cax = ax2c, aspect=90,pad=.001)
cbar2.set_label(L"$s_2^z$")
cbar2.ax.tick_params(labelsize=fs)


ax3.imshow(cur,aspect="auto",extent=[0,dat[end,1],0,2],cmap="jet")
ax3.set_ylabel(L"$V$")
ax3.set_xticklabels([])
ax3.annotate("(c)",(0.013,0.85),xycoords="axes fraction",fontsize=fs)
cmap = plt.cm.ScalarMappable(cmap=cp,norm=plt.Normalize(vmin = minimum(cur), vmax = maximum(cur)))
cbar3 = fig.colorbar(cmap, ax=ax3, cax = ax3c, aspect=90,pad=.001)
cbar3.set_label(L"$I_L - I_R\, [\Gamma_0]$")
cbar3.ax.tick_params(labelsize=fs)

ax4.imshow(X,aspect="auto",extent=[0,dat[end,1],0,2],cmap="jet")
ax4.set_xticklabels([])
ax4.annotate("(d)",(0.015,0.85),xycoords="axes fraction",fontsize=fs)
ax4.set_ylabel(L"$V$")
cmap = plt.cm.ScalarMappable(cmap=cp,norm=plt.Normalize(vmin = minimum(X), vmax = maximum(X)))
cbar4 = fig.colorbar(cmap, ax=ax4, cax = ax4c, aspect=90,pad=.001)
cbar4.set_label(L"$X$")
cbar4.ax.tick_params(labelsize=fs)

ax5.imshow(P,aspect="auto",extent=[0,dat[end,1],0,2],cmap="jet")
ax5.set_xlabel(L"$t\,[t_0^{-1}]$")
ax5.annotate("(e)",(0.015,0.85),xycoords="axes fraction",fontsize=fs)
ax5.set_ylabel(L"$V$")
cmap = plt.cm.ScalarMappable(cmap=cp,norm=plt.Normalize(vmin = minimum(P), vmax = maximum(P)))
cbar5 = fig.colorbar(cmap, ax=ax5, cax = ax5c, aspect=90,pad=.001)
cbar5.set_label(L"$P$")
cbar5.ax.tick_params(labelsize=fs)


plt.savefig("Pops_currents_phonons_imshow_Peierls_1x2_VVar.pdf")





