using PyPlot
using DelimitedFiles
using FnMatch

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 19
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.14
py = 0.58
sx = 0.82
sy = 0.40
dy = 0.05
dx = 0.01

ax1=plt.axes([px,py,sx,sy])
ax2=plt.axes([px,py-dy-sy,sx,sy])
sufP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"
sufM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"

#pre = ["L3/L3","L1x5/L5","L2x3/L6","L3x6/L18"]
#lab = [L"$N=1\times 3$", L"$N=1\times 5$",L"$N=2\times 3$",L"$N=3\times 6$"]
#col = ["k","red","blue","green"]

pre = ["L1x5/L5","L3x6/L18"]
lab = [L"$N=1\times 5$",L"$N=3\times 6$"]
col = ["k","red"]

for (p,c,l) in zip(pre,col,lab)
	datP = readdlm("$(p)$(sufP)")[2:end,:]
	datM = readdlm("$(p)$(sufM)")[2:end,:]
	lP = length(datP[:,1])
	lM = length(datM[:,1])
	lD = min(lP,lM)
	ti = datP[1:lD,1]*fac_t
	IP = (datP[1:lD,2]-datP[1:lD,3])*fac_c/2
	IM = (datM[1:lD,2]-datM[1:lD,3])*fac_c/2
	ax1.plot(ti,IP,color=c,label=l,linewidth=lw-1)
	ax1.plot(ti,IM,color=c,linestyle="dashed",linewidth=lw)
	
	Is = IP+IM
	Id = IP-IM
	SP = abs.(100 * Id ./ Is)
	ax2.plot(ti,SP,color=c,label=l,linewidth=lw-1)
end

dat0P = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
dat0M = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]

l0P = length(dat0P[:,1])
l0M = length(dat0M[:,1])

lD = min(l0P,l0M)
dat0P = dat0P[1:lD,:]
dat0M = dat0M[1:lD,:]

ti = dat0P[:,1]*fac_t
I0P = (dat0P[:,2]-dat0P[:,3])*fac_c/2
I0M = (dat0M[:,2]-dat0M[:,3])*fac_c/2

ax1.plot(ti,I0P,color="pink",label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$",linewidth=lw-1)
ax1.plot(ti,I0M,color="pink",linestyle="dashed",linewidth=lw)

Is = I0P+I0M
Id = I0P-I0M
SP = abs.(100 * Id ./ Is)
ax2.plot(ti,SP,color="pink",linewidth=lw-1)

ax1.set_ylabel("current "*L"$I\,[\mu\mathrm{A}]$")	
ax1.set_xticklabels([])
ax1.set_xlim(0,50)
ax1.legend(loc="lower left")

ax2.set_ylabel("spin selectivity "*L"S\,[\%]")
ax2.set_xlabel("time "*L"$t\,[\mathrm{ps}]$")
ax2.set_xlim(0,50)
ax2.set_yscale("log")
ax2.set_ylim(1e-7,1e2)
plt.savefig("Cur_sp_singleTrag_X0_P0_Nvar.pdf")
	
	
	
	






