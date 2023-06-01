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

px = 0.065
py = 0.09
sx = 0.26
sy = 0.88
dy = 0.03
dx = 0.02

fig = plt.figure(figsize=(9,6))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
cax3 = plt.axes([px+3*sx+2*dx+0.01,py,0.013,sy])

V = 55

ws = 0.002
T0 = 0.04/ws

fac_c_v2 = 1e3 * 0.2433946
fac_t_v2 = 1/(4.135668*1000)
fac_c_v3 = 0.002 * 1e3 * 0.2433946
fac_t_v3 = 1/(ws * 4.135668*1000)

Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\mu\mathrm{A}]$",L"$s_1^z$"]

dat24_v2 = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L2x4/V55/electrons_avg_n1515_comb.dat")
dat24_v3 = readdlm("L2x4/V55/MFHEOM_v3/P/electrons_avg_n41.dat")
dat24_B = readdlm("bundles/V55/L2x4/L8_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V220.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_electrons.dat")


cp = plt.cm.inferno
cmap = plt.cm.ScalarMappable(cmap=cp)
sp = dat24_B[2:end,2:2:end]-dat24_B[2:end,3:2:end]
vx = maximum(sp)
vm = minimum(sp)

ax1.imshow(dat24_v2[2:end,2:2:end]-dat24_v2[2:end,3:2:end],aspect="auto",extent=[1,8,dat24_v2[end,1]*fac_t_v2,dat24_v2[2,1]*fac_t_v2],vmin=vm,vmax=vx,cmap="inferno")
ax1.set_ylim(30,0)
ax1.set_ylabel(L"$t\,[\mathrm{ps}]$")
ax1.set_xlabel("site")
ax2.imshow(dat24_v3[2:end,2:2:end]-dat24_v3[2:end,3:2:end],aspect="auto",extent=[1,8,dat24_v3[end,1]*fac_t_v3,dat24_v3[2,1]*fac_t_v3],vmin=vm,vmax=vx,cmap="inferno")
ax2.set_ylim(30,0)
ax2.set_yticklabels([])
ax2.set_xlabel("site")
im = ax3.imshow(dat24_B[2:end,2:2:end]-dat24_B[2:end,3:2:end],aspect="auto",extent=[1,8,dat24_B[end,1]*fac_t_v3,dat24_B[2,1]*fac_t_v3],vmin=vm,vmax=vx,cmap="inferno")
ax3.set_ylim(30,0)
ax3.set_yticklabels([])
ax3.set_xlabel("site")
ax1.annotate("v1",(0.8,0.07),xycoords="axes fraction",fontsize=fs,color="white")
ax2.annotate("v2",(0.8,0.07),xycoords="axes fraction",fontsize=fs,color="white")
ax3.annotate("bundles",(0.6,0.07),xycoords="axes fraction",fontsize=fs,color="black")

cbar = fig.colorbar(im,ax=ax3,cax = cax3,aspect=90,pad=.001)
cbar.set_label(L"$s_j^z$",labelpad=-8)
cbar.ax.tick_params(labelsize=fs)

plt.savefig("Comparison_MFHEOM_v3_v2_bundles_spinPol_L2x4_V55.pdf")

