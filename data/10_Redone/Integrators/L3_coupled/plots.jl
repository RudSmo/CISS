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
fac_t = 1/(0.002 * 4.135668)

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

datRK = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_RK4_sequential_RK4_VNL_dt0.0001_1_convergence_info.dat") 
datSE1 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_SymplecticEulerStep_RK4_VNL_dt0.0001_1_convergence_info.dat")
datSE2 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_SymplecticEulerStep_RK4_VNL_dt1.0e-5_1_convergence_info.dat") 
datSV  = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_SymplecticVerletStep_RK4_VNL_dt0.0001_1_convergence_info.dat")
datSR3 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_SymplecticRuth3Step_RK4_VNL_dt0.0001_1_convergence_info.dat")
datSR4 = readdlm("L3_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_SymplecticRuth4Step_RK4_VNL_dt0.0001_1_convergence_info.dat")

axes = [ax1,ax2,ax3]
col = ["blue","green","red","turquoise","purple","black"]

ax1.plot(datRK[2:end,1]*fac_t,abs.(datRK[2:end,3] .- datRK[2,3]) .+ 1e-9,label="RK4",linewidth=lw,color=col[1])
ax1.plot(datSE1[2:end,1]*fac_t,abs.(datSE1[2:end,3] .- datSE1[2,3]) .+ 1e-9,label="symp. Euler, "*L"$dt=0.0001$",linewidth=lw,color=col[2])
ax1.plot(datSE2[2:end,1]*fac_t,abs.(datSE2[2:end,3] .- datSE2[2,3]) .+ 1e-9,label="symp. Euler, "*L"$dt=0.00001$",linewidth=lw,color=col[3])
ax1.plot(datSV[2:end,1]*fac_t,abs.(datSV[2:end,3] .- datSV[2,3]) .+ 1e-9,label="symp. Verlet",linewidth=lw,linestyle="dashed",color=col[4])
ax1.plot(datSR3[2:end,1]*fac_t,abs.(datSR3[2:end,3] .- datSR3[2:end,3]) .+ 1e-9,label="symp. Ruth",linewidth=lw,color=col[5])
ax1.plot(datSR4[2:end,1]*fac_t,abs.(datSR4[2:end,3] .- datSR4[2:end,3]) .+ 1e-9,label="symp. Forest-Ruth",linewidth=lw,linestyle="dashed",color=col[6])
ax1.set_yscale("log")
ax1.set_xlim(0,1200)
ax1.set_ylim(1e-10,1e-2)
ax1.set_xticklabels([])
annx = 0.02
anny = 0.9
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel(L"$\vert E_\mathrm{tot}(t)-E_\mathrm{tot}(t=0)\vert+10^{-9}$")


ind = [2,4]
ax2.plot(datRK[2:end,1]*fac_t,datRK[2:end,2],linewidth=lw,color=col[1],label="RK4")
ax2.plot(datSE1[2:end,1]*fac_t,datSE1[2:end,2],linewidth=lw,color=col[2],label="symp. Euler, "*L"$dt=0.0001$")
ax2.plot(datSE2[2:end,1]*fac_t,datSE2[2:end,2],linewidth=lw,color=col[3],label="symp. Euler, "*L"$dt=0.00001$")
ax2.plot(datSV[2:end,1]*fac_t,datSV[2:end,2],linewidth=lw,linestyle="dashed",color=col[4],label="symp. Verlet")
ax2.plot(datSR3[2:end,1]*fac_t,datSR3[2:end,2],linewidth=lw,color=col[5],label="symp. Ruth")
ax2.plot(datSR4[2:end,1]*fac_t,datSR4[2:end,2],linewidth=lw,linestyle="dashed",color=col[6],label="symp. Forest-Ruth")	
ax2.set_xlim(0,1200)
ax2.set_ylabel(L"$\mathrm{tr}\{\rho(t)\}$")
ax2.set_xticklabels([])
ax2.annotate("(b)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax2.legend(loc="center right",framealpha=0.8)

fac_4 = 1e16
ax3.plot(datRK[2:end,1]*fac_t,datRK[2:end,4]*fac_4,linewidth=lw,color=col[1])
ax3.plot(datSE1[2:end,1]*fac_t,datSE1[2:end,4]*fac_4,linewidth=lw,color=col[2])
ax3.plot(datSE2[2:end,1]*fac_t,datSE2[2:end,4]*fac_4,linewidth=lw,color=col[3])
ax3.plot(datSV[2:end,1]*fac_t,datSV[2:end,4]*fac_4,linewidth=lw,linestyle="dashed",color=col[4])
ax3.plot(datSR3[2:end,1]*fac_t,datSR3[2:end,4]*fac_4,linewidth=lw,color=col[5])
ax3.plot(datSR4[2:end,1]*fac_t,datSR4[2:end,4]*fac_4,linewidth=lw,linestyle="dashed",color=col[6])
ax3.set_xlim(0,1200)
ax3.set_ylabel(L"$\dot{W}_\mathrm{el-vib}+\dot{E}_\mathrm{vib} \, [\times 10^{-16}]$")
ax3.annotate("(c)",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")

plt.savefig("L3_coupled_iso_integrators.pdf")

