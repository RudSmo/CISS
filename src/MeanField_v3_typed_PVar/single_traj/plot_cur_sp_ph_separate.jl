using PyPlot
using DelimitedFiles
using FnMatch

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 21
lw = 4
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
rc("legend",fontsize=fs-5)
px = 0.16
py = 0.18
sx = 0.81
sy = 0.78
dy = 0.04
dx = 0.01

fig=plt.figure(figsize=(6,4))
ax=plt.axes([px,py,sx,sy])
c = "tab:blue"
l = ""
lws = lw-2

datP=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
datM=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]


lP = length(datP[:,1])
lM = length(datM[:,1])
lD = min(lP,lM)
ti = datP[1:lD,1]*fac_t
IP = (datP[1:lD,2]-datP[1:lD,3])*fac_c/2
IM = (datM[1:lD,2]-datM[1:lD,3])*fac_c/2
ax.plot(ti,IP,color=c,label=l,linewidth=lws)
ax.plot(ti,IM,color=c,linestyle="dashed",linewidth=lws)
ax.set_ylabel("current "*L"$I\,[\mu\mathrm{A}]$",labelpad=15)	
ax.set_xlabel("time "*L"$t\, [\mathrm{ps}]$")
ax.set_xlim(0,50)
plt.savefig("L5_init_X0_P0_current_PM.pdf")

	
plt.clf()
fig=plt.figure(figsize=(6,4))
ax=plt.axes([px,py,sx,sy])
Is = IP+IM
Id = IP-IM
SP = 100 * Id ./ Is
ax.plot(ti,SP,color=c,label=l,linewidth=lws)
ax.set_ylabel("spin selectivity "*L"S\,[\%]")
ax.set_xlim(0,50)
ax.set_xlabel("time "*L"$t\, [\mathrm{ps}]$")
ax.set_ylim(-80,100)
plt.savefig("L5_init_X0_P0_spin_selectivity.pdf")


plt.clf()
fig=plt.figure(figsize=(6,4))
ax=plt.axes([px,py,sx,sy])
datXP=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
datXM=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
lXP = length(datXP[:,1])
lXM = length(datXM[:,1])
lf = min(lXP,lXM)
datXP = datXP[1:lf,:]
datXM = datXM[1:lf,:]
tim = datXP[:,1]*fac_t
XP = datXP[:,2]+datXP[:,4]+datXP[:,6]+datXP[:,8]
XM = datXM[:,2]+datXM[:,4]+datXM[:,6]+datXM[:,8]
dX = XP-XM
ax.plot(tim,dX,color=c,linewidth=lw-1,label=l)
ax.set_xlim(0,50)
ax.set_ylim(-10,10)
ax.set_ylabel("displacement polar. "*L"$S_X$")
ax.set_xlabel("time "*L"$t\,[\mathrm{ps}]$")
plt.savefig("L5_init_X0_P0_Xdif.pdf")
	
	
	
	






