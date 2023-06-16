using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002 * 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002 * 4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,6))
ax = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+0.1,py+0.4,sx-0.45,sy-0.48])

num_sites = 3
num_laps = 2
N = num_laps * num_sites
P = 0.5      


dat = readdlm("data0P/L6_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_currents_avg.dat")
ax.plot(dat[2:end,1]*fac_t*1e-3,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw,label="Trajectory "*L"$n=25, P = 0.5$",color="red")
ax2.plot(dat[2:end,1]*fac_t*1e-3,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw,label="Trajectory "*L"$n=25, P = 0.5$",color="red")

dat = readdlm("data0M/L6_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_currents_avg.dat")
ax.plot(dat[2:end,1]*fac_t*1e-3,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw,label="Trajectory "*L"$n=25, P = -0.5$",color="blue")
ax2.plot(dat[2:end,1]*fac_t*1e-3,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw,label="Trajectory "*L"$n=25, P = -0.5$",color="blue")

dat = readdlm("currents_avg_n71_P.dat")
ax.plot(dat[:,1]*fac_t*1e-3,(dat[:,2]-dat[:,3])*fac_c,linewidth=lw,label="Average "*L"$n=525, P = 0.5$",color="k")	
ax2.plot(dat[:,1]*fac_t*1e-3,(dat[:,2]-dat[:,3])*fac_c,linewidth=lw,label="Average "*L"$n=525, P = 0.5$",color="k")

dat = readdlm("currents_avg_n71_M.dat")
ax.plot(dat[:,1]*fac_t*1e-3,(dat[:,2]-dat[:,3])*fac_c,linewidth=lw,label="Average "*L"$n=525, P = -0.5$",color="orange")
ax2.plot(dat[:,1]*fac_t*1e-3,(dat[:,2]-dat[:,3])*fac_c,linewidth=lw,label="Average "*L"$n=525, P = -0.5$",color="orange")
		
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
#ax.set_xlim(-0.1,30.1)
ax.legend()
ax2.set_xlim(82,84)
ax2.set_ylim(0.4,0.8)
ax.annotate(L"$L=2\times 3,\ V=0.48\, \mathrm{V}$",(0.7,0.55),xycoords="axes fraction",fontsize=fs)
plt.savefig("L$(num_laps)x$(num_sites)_P$(P)_SampleVar_current.pdf")
