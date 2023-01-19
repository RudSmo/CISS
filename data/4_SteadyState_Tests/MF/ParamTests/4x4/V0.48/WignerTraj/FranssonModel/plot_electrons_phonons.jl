using LinearAlgebra
using DelimitedFiles
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
lw = 3
al = 0.4
fac_t = 1/(4.135668*1e3)
fac_c = 1e3 * 2.8494561
fac_ph = 0.002

#plot = "electrons"
plot = "spin"
#phonons = false
phonons = true

N = 16
V = 0.48
tPh = 0.004
#tPh = 0.0
tso = 0.001
tsoPh = 0.0001
#tsoPh = 0.0
Tend = 37500
#Tend = 6250
ver = 2
P = 0.0

dat_el = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P$(P)_X0_-10_dt0.05_$(ver)_electrons.dat")
dat_ph = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P$(P)_X0_-10_dt0.05_$(ver)_phonons.dat")
dat_cur = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P$(P)_X0_-10_dt0.05_$(ver)_currents.dat")

fig = plt.figure(figsize=(10,8))
px = 0.08
py = 0.07
sx1 = 0.19
sx2 = 0.49
sy = 0.83
dx = 0.02
dy = 0.01

ax1 = plt.axes([px,py,sx1,sy])
ax2 = plt.axes([px+sx1+dx,py,sx2,sy])
ax2_c = plt.axes([px+sx1+dx,py+sy+dy,sx2,0.02])

ax3 = plt.axes([px+sx1+dx+sx2+dx,py,sx1,sy])

ax1.plot(dat_cur[2:end,2]*fac_c,dat_cur[2:end,1]*fac_t,color="k",linewidth=lw)
ax1.set_ylim(0,Tend*fac_t)
ax1.set_ylabel(L"$t\,[\mathrm{ps}]$")
ax1.set_xlabel(L"$I_L\,[\mu\mathrm{A}]$")
ax1.set_xticks([-5,0,5,10])
ax1.set_xlim(-5,10)

if plot == "electrons"
	Nd = dat_el[2:end,2:2:end] + dat_el[2:end,3:2:end]
else
	Nd = (dat_el[2:end,2:2:end] - dat_el[2:end,3:2:end]) * 10^3
end

cp = plt.cm.jet
vmn = minimum(Nd)
vmx = maximum(Nd)
norm = plt.matplotlib.colors.Normalize(vmin=vmn, vmax=vmx)
cmap = plt.cm.ScalarMappable(cmap=cp,norm=norm)

ax2.imshow(Nd,aspect="auto",extent=[1,N,Tend*fac_t,0],vmin=vmn,vmax=vmx,cmap="jet")

if phonons 
	for i in 1:N
	
		pos = i * ones(length(dat_ph[2:end,1]))	
		if i == 1
			ax2.plot(pos - fac_ph * dat_ph[2:end,2*i], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw,alpha=al)
		elseif i == N
			ax2.plot(pos + fac_ph * dat_ph[2:end,2*(i-1)], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw,alpha=al)
		else
			ax2.plot(pos - fac_ph * dat_ph[2:end,2*(i-1)] + fac_ph * dat_ph[2:end,2*i], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw,alpha=al)
		end
	
	
	end
end

ax2.set_ylim(0,Tend*fac_t)
ax2.set_xlim(0.5,N+0.5)
#ax2.set_yticks(Array(0.0:0.5:3.0))
ax2.set_yticklabels([])
ax2.set_xlabel("site")
cbar = fig.colorbar(cmap, ax=ax2, cax = ax2_c,pad=.001,aspect=90,orientation="horizontal")
if plot == "electrons"
	cbar.set_label(L"$n$",labelpad=-47)
else
	cbar.set_label(L"$s^z\,[\times 10^{-3}]$",labelpad=-49)
end
cbar.ax.tick_params(labelsize=fs,labelbottom=false,labeltop=true)

ax3.set_ylim(0,Tend*fac_t)
ax3.set_xlim(-15,1)
ax3.set_xticks([-15,-10,-5,0])
ax3.set_yticklabels([])
ax3.plot(dat_cur[2:end,3]*fac_c,dat_cur[2:end,1]*fac_t,color="k",linewidth=lw)
ax3.set_xlabel(L"$I_R\,[\mu\mathrm{A}]$")

plt.savefig("L$(N)_V$(V)_tPh$(tPh)_tso$(tso)_tsoPh$(tsoPh)_P$(P)_v$(ver)_electrons_currents_phonons.pdf")


