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
dat11 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_currents.dat") 
dat12 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticEulerStep_RK4_VNL!_dt0.0001_1_currents.dat")
dat13 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_RK4_sequential_RK4_VNL!_dt0.0001_1_currents.dat") 

# electrons 
dat21 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_electrons.dat")
dat22 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticEulerStep_RK4_VNL!_dt0.0001_1_electrons.dat")
dat23 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_RK4_sequential_RK4_VNL!_dt0.0001_1_electrons.dat")

# phonons
dat31 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0001_1_phonons.dat")
dat32 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticEulerStep_RK4_VNL!_dt0.0001_1_phonons.dat")
dat33 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_RK4_sequential_RK4_VNL!_dt0.0001_1_phonons.dat")



time = dat11[2:end,1]
dat11 = dat11[2:end,2:end]
dat12 = dat12[2:end,2:end]
dat13 = dat13[2:end,2:end]
dat21 = dat21[2:end,2:end]
dat22 = dat22[2:end,2:end]
dat23 = dat23[2:end,2:end]
dat31 = dat31[2:end,2:end]
dat32 = dat32[2:end,2:end]
dat33 = dat33[2:end,2:end]

axes = [ax1,ax2,ax3]
datas = [[dat11,dat12,dat13],[dat21,dat22,dat23],[dat31,dat32,dat33]]

labels = ["Symplectic Ruth 4th","Symplectic Euler 1st","Runge-Kutta 4th"]
colors = ["blue","green","red"]

ylabels = [L"$I\,[\Gamma_0]$",L"$s_1^z$",L"$X_1$"]
annotations = ["(a)","(b)","(c)"]
linewidths = [lw,lw,lw]
ws = 0.002
T0 = 0.04/ws
fac_t = 1/(ws * 4.135668)
Gamma = T0/4
annx = 0.02
anny = 0.9


for (ax,data,ylab,ann) in zip(axes,datas,ylabels,annotations)
	for (col,dat,lab,lw) in zip(colors,data,labels,linewidths)
		if ax == axes[1]
			ax.plot(time*fac_t,(dat[:,1]-dat[:,2])/Gamma,label=lab,linewidth=lw,color=col)
		elseif ax == axes[2]
			ax.plot(time*fac_t,dat[:,1]-dat[:,2],label=lab,linewidth=lw,color=col)
		else
			ax.plot(time*fac_t,dat[:,1],label=lab,linewidth=lw,color=col)
		end
	end
	ax.set_ylabel(ylab,labelpad=3)
	ax.set_xlim(0,6000)
	ax.annotate(ann,(annx,anny),xycoords="axes fraction",fontsize=fs)
end

ax1.set_ylim(0.25,0.6)
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel(L"$I \, [\Gamma]$")
ax1.legend(loc="center right",framealpha=0.8)
ax2.set_ylim(0.105,0.14)
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax3.set_ylim(-12,9)
plt.savefig("L3_convergence_integrators.pdf")

