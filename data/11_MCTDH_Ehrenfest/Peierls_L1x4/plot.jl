using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.69
sx = 0.87
sy = 0.27
dy = 0.03

fig = plt.figure(figsize=(9,9))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx,sy])

ws = 0.001

fac_t = 1/(ws*4.135668)

annx = 0.02
anny = 0.9

datM = readdlm("MCTDH/MCTDH_L4_Peierls.dat")
datE = readdlm("Ehrenfest/undamped/L4_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V300.0_w00.001_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_electrons_avg.dat")
datEd = readdlm("Ehrenfest/damped/damp0.5/L4_tPh_2.0_thPd_2.0_SOC0.0_SOCPh0.0_V300.0_w00.001_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_electrons_avg.dat")

col = ["purple","blue","green","red"]
for i in 2:5
	ax1.plot(datM[2:end,1],datM[2:end,i],linewidth=lw,color=col[i-1],label="site $(i-1)")
end

ax1.legend()
ax1.set_xticklabels([])

ax1.annotate("ML-MCTDH",(0.8,0.07),xycoords="axes fraction",fontsize=fs)
ax1.set_ylabel("Population")
annx = 0.013
anny = 0.88
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",fontsize=fs)
for i in 2:2:8
        ax2.plot(datE[2:end,1]*fac_t,datE[2:end,i],linewidth=lw,color=col[Int(i/2)],label="site $(i-1)")
end

ax2.annotate("Ehrenfest "*L"$(n_\mathrm{traj}=500)$",(0.7,0.85),xycoords="axes fraction",fontsize=fs)
ax2.set_ylabel("Population")
ax2.set_xticklabels([])
ax2.annotate("(b)",(annx,anny),xycoords="axes fraction",fontsize=fs)
for i in 2:2:8
        ax3.plot(datEd[2:end,1]*fac_t,datEd[2:end,i],linewidth=lw,color=col[Int(i/2)],label="site $(i-1)")
end

ax3.annotate("Ehrenfest "*L"$(n_{\mathrm{traj}}=500)$"*" with damping "*L"$\gamma = 0.002$",(0.4,0.85),xycoords="axes fraction",fontsize=fs)
ax3.set_ylabel("Population")
ax3.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax3.annotate("(c)",(annx,anny),xycoords="axes fraction",fontsize=fs)


plt.savefig("Peierls_L1x4_MCTDH_v_Ehrenfest.pdf")

