using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 1.5
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.7
sx = 0.87
sy = 0.28
dy = 0.03

fig = plt.figure(figsize=(9,9))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
V = 140.0
axes = [ax1,ax2,ax3]
colors = ["blue","green","red","orange","black","purple"]
annotations = ["(a)","(b)","(c)"]
ws = 0.002
T0 = 0.04/ws
fac_t = 1/(ws * 4.135668 * 1000)
Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\Gamma_0]$",L"$s_1^z$",L"$X$"]

for (i,col) in zip(0:5,colors)
	dat1 = readdlm("data$(i)/L2_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V$(V)_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_currents_avg.dat") 
	dat2 = readdlm("data$(i)/L2_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V$(V)_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_electrons_avg.dat")
	dat3 = readdlm("data$(i)/L2_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V$(V)_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_phonons_avg.dat")
	datas = [dat1,dat2,dat3]
	for (ax,dat,ylab,ann) in zip(axes,datas,ylabels,annotations)
		if ax == axes[1]
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2]/Gamma,linewidth=lw,color=col)
		elseif ax == axes[2]
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2]-dat[2:end,3],linewidth=lw,color=col)
		else
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2],linewidth=lw,color=col)
		end
		ax.set_ylabel(ylab)
		ax.set_xlim(0,180)
		ax.annotate(ann,(annx,anny),xycoords="axes fraction",fontsize=fs)
	end
end

ax1.set_xticklabels([])
ax1.set_ylim(0.0,0.05)
ax2.set_xticklabels([])
ax2.set_ylim(0.0,0.02)
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel(L"$I \, [\Gamma]$")
#ax1.legend(loc="lower right",framealpha=0.8)

ax3.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax3.set_ylim(-7,7)

plt.savefig("L2_dynamics_V$(V).pdf")

