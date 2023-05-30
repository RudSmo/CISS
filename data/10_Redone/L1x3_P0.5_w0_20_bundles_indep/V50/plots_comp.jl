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
ax11 = plt.axes([0.23,0.79,0.42,0.15])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax21 = plt.axes([0.25,0.43,0.6,0.15])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])
V = 50
axes = [ax1,ax2,ax3]
ann = ["(a)","(b)","(c)"]
ws = 0.02
T0 = 0.04/ws
fac_t = 1/(ws * 4.135668 * 1000)
Gamma = T0/4
annx = 0.02
anny = 0.9
ylab = [L"$I\,[\Gamma_0]$",L"$s_1^z$",L"$X$"]

dat1 = readdlm("n500000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_currents.dat")
dat2 = readdlm("n500000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_electrons.dat")
dat3 = readdlm("n500000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_phonons.dat")

dat11 = readdlm("n10000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_currents.dat")
dat21 = readdlm("n10000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_electrons.dat")
dat31 = readdlm("n10000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_phonons.dat")

dat12 = readdlm("n100000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_currents.dat")
dat22 = readdlm("n100000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_electrons.dat")
dat32 = readdlm("n100000_bundles/L3_tPh_0.2_thPd_0.2_SOC0.05_SOCPh0.005000000000000001_V20.0_w00.02_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_phonons.dat")

dati1 = readdlm("n15150_indep/currents_avg_n101.dat")
dati2 = readdlm("n15150_indep/electrons_avg_n101.dat")
dati3 = readdlm("n15150_indep/phonons_avg_n101.dat")

ax1.plot(dati1[2:end,1]*fac_t,(dati1[2:end,2]-dati1[2:end,3])/Gamma,linewidth=lw,color="blue",label="independent "*L"$(n=15150)$")
ax1.plot(dat11[2:end,1]*fac_t,(dat11[2:end,2]-dat11[2:end,3])/Gamma,linewidth=lw,color="green",label="bundles "*L"$(n=10000)$")
ax1.plot(dat12[2:end,1]*fac_t,(dat12[2:end,2]-dat12[2:end,3])/Gamma,linewidth=lw,color="black",label="bundles "*L"$(n=100000)$",alpha=0.8)
ax1.plot(dat1[2:end,1]*fac_t,(dat1[2:end,2]-dat1[2:end,3])/Gamma,linewidth=lw,color="red",label="bundles "*L"$(n=500000)$",alpha=0.6)

ax1.set_ylabel(ylab[1])
ax1.annotate(ann[1],(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_xticklabels([])
ax1.set_xlim(0,1.7)
ax1.set_ylim(0.22,0.3)
ax1.legend(loc="best")

ax11.plot(dati1[2:end,1]*fac_t,(dati1[2:end,2]-dati1[2:end,3])/Gamma,linewidth=lw,color="blue",label="independent "*L"$(n=15150)$")
ax11.plot(dat11[2:end,1]*fac_t,(dat11[2:end,2]-dat11[2:end,3])/Gamma,linewidth=lw,color="green",label="bundles "*L"$(n=10000)$")
ax11.plot(dat12[2:end,1]*fac_t,(dat12[2:end,2]-dat12[2:end,3])/Gamma,linewidth=lw,color="black",label="bundles "*L"$(n=100000)$",alpha=0.8)
ax11.plot(dat1[2:end,1]*fac_t,(dat1[2:end,2]-dat1[2:end,3])/Gamma,linewidth=lw,color="red",label="bundles "*L"$(n=500000)$",alpha=0.6)
ax11.set_xlim(0.4,1.7)
ax11.set_ylim(0.233,0.235)

ax2.plot(dati2[2:end,1]*fac_t,dati2[2:end,2]-dati2[2:end,3],linewidth=lw,color="blue")
ax2.plot(dat21[2:end,1]*fac_t,dat21[2:end,2]-dat21[2:end,3],linewidth=lw,color="green",label="bundles "*L"$(n=10000)$")
ax2.plot(dat22[2:end,1]*fac_t,dat22[2:end,2]-dat22[2:end,3],linewidth=lw,color="black",label="bundles "*L"$(n=100000)$",alpha=0.8)
ax2.plot(dat2[2:end,1]*fac_t,dat2[2:end,2]-dat2[2:end,3],linewidth=lw,color="red",alpha=0.6)
ax2.set_ylabel(ylab[2])
ax2.annotate(ann[2],(annx,anny),xycoords="axes fraction",fontsize=fs)
ax2.set_xticklabels([])
#ax2.set_ylim(0.0,0.02)
ax2.set_xlim(0,1.7)

ax21.plot(dati2[2:end,1]*fac_t,dati2[2:end,2]-dati2[2:end,3],linewidth=lw,color="blue")
ax21.plot(dat21[2:end,1]*fac_t,dat21[2:end,2]-dat21[2:end,3],linewidth=lw,color="green",label="bundles "*L"$(n=10000)$")
ax21.plot(dat22[2:end,1]*fac_t,dat22[2:end,2]-dat22[2:end,3],linewidth=lw,color="black",label="bundles "*L"$(n=100000)$",alpha=0.8)
ax21.plot(dat2[2:end,1]*fac_t,dat2[2:end,2]-dat2[2:end,3],linewidth=lw,color="red",alpha=0.6)
ax21.set_xlim(0.4,1.7)
ax21.set_ylim(0.065,0.0657)

ax3.plot(dati3[2:end,1]*fac_t,dati3[2:end,2],linewidth=lw,color="blue")
ax3.plot(dat31[2:end,1]*fac_t,dat31[2:end,2],linewidth=lw,color="green",label="bundles "*L"$(n=10000)$")
ax3.plot(dat32[2:end,1]*fac_t,dat32[2:end,2],linewidth=lw,color="black",label="bundles "*L"$(n=100000)$",alpha=0.8)
ax3.plot(dat3[2:end,1]*fac_t,dat3[2:end,2],linewidth=lw,color="red",alpha=0.6)
ax3.set_ylabel(ylab[3])
ax3.annotate(ann[3],(annx,anny),xycoords="axes fraction",fontsize=fs)
ax3.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax3.set_ylim(-0.13,0.05)
ax3.set_xlim(0,1.7)

plt.savefig("L1x3_w0_20_comparison_bundles_indep_V$(V).pdf")

