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
fac_t = 1/(4.135668*1e3)
fac_c = 1e3 * 2.8494561
fac_ph = 0.001

N = 16
V = 0.48
tPh = 0.004
tso = 0.0
tsoPh = 0.0
Tend = 12500


dat_el = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_electrons.dat")
dat_ph = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_phonons.dat")
dat_cur = readdlm("L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat")

fig = plt.figure(figsize=(10,8))
px = 0.08
py = 0.08
sx1 = 0.19
sx2 = 0.49
sy = 0.9
dx = 0.02

ax1 = plt.axes([px,py,sx1,sy])
ax2 = plt.axes([px+sx1+dx,py,sx2,sy])
ax3 = plt.axes([px+sx1+dx+sx2+dx,py,sx1,sy])


ax1.plot(dat_cur[2:end,2]*fac_c,dat_cur[2:end,1]*fac_t,color="k",linewidth=lw)
ax1.set_ylim(0,Tend*fac_t)
ax1.set_ylabel(L"$t\,[\mathrm{ps}]$")
ax1.set_xlabel(L"$I_L\,[\mathrm{mA}]$")
ax1.set_xticks([-5,0,5,10])
ax1.set_xlim(-5,10)

ax2.imshow(dat_el[2:end,2:2:end] .+ dat_el[2:end,3:2:end],aspect="auto",extent=[1,N,Tend*fac_t,0],origin="lower")


#ax2.plot(ones(length(dat_ph[2:end,1])),dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)


for i in 1:N
#	ax2.plot((i+1)*ones(length(dat_ph[2:end,1])) + 2*0.003*dat_ph[2:end,2*i],dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
	pos = i * ones(length(dat_ph[2:end,1]))	
	if i == 1
		ax2.plot(pos - fac_ph * dat_ph[2:end,2*i], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
	elseif i == N
		ax2.plot(pos + fac_ph * dat_ph[2:end,2*(i-1)], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
	else
		ax2.plot(pos - fac_ph * dat_ph[2:end,2*(i-1)] + fac_ph * dat_ph[2:end,2*i], dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
	end


#	ax2.plot(dat_ph[2:end,2*i]-2*dat_ph[2:end,2*i]*0.003/2 .+ i,dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
#	ax2.plot(N - 2*dat_ph[2:end,2*i]*0.003/2 .+ i,dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
#	ax2.plot(N - 2*dat_ph[2:end,2*i]*0.003/2 .+ i,dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
#	ax2.plot(dat_ph[2:end,2*(i+1)]*0.003/2 .+ i+1,dat_ph[2:end,1]*fac_t,color="k",linewidth=lw)
end

ax2.set_ylim(0,Tend*fac_t)
ax2.set_xlim(0.5,N+0.5)
ax2.set_yticks(Array(0.0:0.5:3.0))
ax2.set_yticklabels([])
ax2.set_xlabel("site")

ax3.set_ylim(0,Tend*fac_t)
ax3.set_xlim(-15,1)
ax3.set_xticks([-15,-10,-5,0])
ax3.set_yticklabels([])
ax3.plot(dat_cur[2:end,3]*fac_c,dat_cur[2:end,1]*fac_t,color="k",linewidth=lw)
ax3.set_xlabel(L"$I_R\,[\mathrm{mA}]$")

plt.savefig("L$(N)_V$(V)_tPh$(tPh)_tso$(tso)_tsoPh$(tsoPh)_electrons_currents_phonons.pdf")


