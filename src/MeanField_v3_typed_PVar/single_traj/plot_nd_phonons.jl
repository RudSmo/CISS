using PyPlot
using DelimitedFiles
using FnMatch

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 13
lw = 2
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
rc("legend",fontsize=fs-5)
px = 0.07
py = 0.08
sx1 = 0.2
sx2= 0.44
sy = 0.82

dy = 0.04
dx = 0.04
fig = plt.figure(figsize=(8,6))

ax1 = plt.axes([px,py,sx1,sy])
ax2 = plt.axes([px+sx1+dx,py,sx2,sy])
ax22 = plt.axes([px+sx1+dx,py+sy+dy,sx2,0.02])
ax3 = plt.axes([px+sx1+sx2+2*dx,py,sx1,sy])

datC = readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")
datX = readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")

ax1.plot(datC[2:end,2],fac_t*datC[2:end,1],color="k",linewidth=lw)

datN = readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
nd = (datN[2:end,2:2:end] .+ datN[2:end,3:2:end] ) ./ 2
cp = plt.cm.plasma
cmap = plt.cm.ScalarMappable(cmap=cp)
mn = minimum(nd)
mx = maximum(nd)
norm = plt.Normalize(mn, mx)

im = ax2.imshow(nd[end:-1:1,:],aspect="auto",extent=[0.5,5.5,0,datC[end,1]*fac_t],vmin=mn,vmax=mx,cmap="inferno",interpolation="none")
cbar = fig.colorbar(im, ax=ax2, cax = ax22, orientation="horizontal", pad=.001)
cbar.set_label("population")
cbar.ax.xaxis.set_label_position("top")
ax2.set_xlabel("site")
#ax2.set_yticks([],[])
ax1.set_ylim(0,50)
ax1.set_xlim(-2,10)
ax2.set_ylim(0,50)

fac_X = 0.008
datX[2:end,2:end] = datX[2:end,2:end]*fac_X
ax2.plot(1 .+ datX[2:end,2]/2,datX[2:end,1]*fac_t,color="k")
ax2.plot(2 .- datX[2:end,2]/2+datX[2:end,4]/2,datX[2:end,1]*fac_t,color="k")
ax2.plot(3 .- datX[2:end,4]/2+datX[2:end,6]/2,datX[2:end,1]*fac_t,color="k")
ax2.plot(4 .- datX[2:end,6]/2+datX[2:end,8]/2,datX[2:end,1]*fac_t,color="k")
ax2.plot(5 .- datX[2:end,8]/2,datX[2:end,1]*fac_t,color="k")

ax2.set_yticklabels([])
ax3.set_yticklabels([])
ax3.set_ylim(0,50)
ax3.set_xlim(-10,2)
ax3.set_xlabel(L"$I_R\, [\mu\mathrm{A}]$")
ax3.plot(datC[2:end,3],fac_t*datC[2:end,1],color="k",linewidth=lw)


ax1.set_ylabel("time "*L"t \,[\mathrm{ps}]")
ax1.set_xlabel(L"$I_L\, [\mu\mathrm{A}]$")



plt.savefig("Cur_nd_phonons_L1x5.pdf")
	
	
	
	






