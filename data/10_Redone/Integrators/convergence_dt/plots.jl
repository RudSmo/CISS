using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
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

# currents 
dat11 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_currents.dat") 
dat12 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_currents.dat")
dat13 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt1.0e-5_1_currents.dat") 

# electrons 
dat21 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_electrons.dat")
dat22 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_electrons.dat")
dat23 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt1.0e-5_1_electrons.dat")

# phonons
dat31 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_phonons.dat")
dat32 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_phonons.dat")
dat33 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V80.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt1.0e-5_1_phonons.dat")

axes = [ax1,ax2,ax3]
datas = [[dat11,dat12,dat13],[dat21,dat22,dat23],[dat31,dat32,dat33]]
labels = [L"$dt=0.024\,\mathrm{fs}$",L"$dt=0.012\,\mathrm{fs}$",L"$dt=0.0012\,\mathrm{fs}$"]
colors = ["blue","green","red"]
ylabels = [L"$I\,[\Gamma_0]$",L"$s_1^z$",L"$X$"]
annotations = ["(a)","(b)","(c)"]
ws = 0.002
T0 = 0.04/ws
fac_t = 1/(ws * 4.135668)
Gamma = T0/4
annx = 0.02
anny = 0.9
for (ax,data,ylab,ann) in zip(axes,datas,ylabels,annotations)
	for (col,dat,lab) in zip(colors,data,labels)
		if ax == axes[1]
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2]/Gamma,label=lab,linewidth=lw,color=col)
		elseif ax == axes[2]
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2]-dat[2:end,3],label=lab,linewidth=lw,color=col)
		else
			ax.plot(dat[2:end,1]*fac_t,dat[2:end,2],label=lab,linewidth=lw,color=col)
		end
	end
	ax.set_ylabel(ylab)
	ax.set_xlim(0,15000)
	ax.annotate(ann,(annx,anny),xycoords="axes fraction",fontsize=fs)
end

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel(L"$I \, [\Gamma]$")
ax1.legend(loc="center right",framealpha=0.8)

ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")

plt.savefig("L3_convergence.pdf")

