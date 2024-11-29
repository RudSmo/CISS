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
px = 0.14
py = 0.705
sx = 0.82
sy = 0.27
dy = 0.04
dx = 0.01

fig=plt.figure(figsize=(8,10))
ax1=plt.axes([px,py,sx,sy])
#ax11=plt.axes([0.21,0.74,0.36,0.16])
ax2=plt.axes([px,py-dy-sy,sx,sy])
ax22=plt.axes([0.25,0.55,0.36,0.1])
ax3=plt.axes([px,py-2*(dy+sy),sx,sy])

sufP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"
sufM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"

#pre = ["L3/L3","L1x5/L5","L2x3/L6","L3x6/L18"]
#lab = [L"$N=1\times 3$", L"$N=1\times 5$",L"$N=2\times 3$",L"$N=3\times 6$"]
#col = ["k","red","blue","green"]

pre = ["L1x5_init_random/L5","L1x5/L5"]
lab = [L"$N=1\times 5$"*" random init",L"$N=1\times 5, X_j(0)=P_j(0)=0$"]
col = ["tab:blue","tab:orange","tab:green","tab:red"]
LW = [lw-1,lw-1,lw-1,lw-2]
for (p,c,l,lws) in zip(pre,col,lab,LW)
	datP = readdlm("$(p)$(sufP)")[2:end,:]
	datM = readdlm("$(p)$(sufM)")[2:end,:]
	lP = length(datP[:,1])
	lM = length(datM[:,1])
	lD = min(lP,lM)
	ti = datP[1:lD,1]*fac_t
	IP = (datP[1:lD,2]-datP[1:lD,3])*fac_c/2
	IM = (datM[1:lD,2]-datM[1:lD,3])*fac_c/2
	ax1.plot(ti,IP,color=c,label=l,linewidth=lws)
	ax1.plot(ti,IM,color=c,linestyle="dashed",linewidth=lws)
#	ax11.plot(ti,IP,color=c,label=l,linewidth=lw-1)
#	ax11.plot(ti,IM,color=c,linestyle="dashed",linewidth=lw-1)
	
	Is = IP+IM
	Id = IP-IM
	SP = 100 * Id ./ Is
	ax2.plot(ti,SP,color=c,label=l,linewidth=lws)
	ax22.plot(ti,SP,color=c,label=l,linewidth=lw-1)
end

#dat0P = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
#dat0M = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
#
#l0P = length(dat0P[:,1])
#l0M = length(dat0M[:,1])
#
#lD = min(l0P,l0M)
#dat0P = dat0P[1:lD,:]
#dat0M = dat0M[1:lD,:]
#
#ti = dat0P[:,1]*fac_t
#I0P = (dat0P[:,2]-dat0P[:,3])*fac_c/2
#I0M = (dat0M[:,2]-dat0M[:,3])*fac_c/2
#
#ax1.plot(ti,I0P,color=col[3],label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$",linewidth=lw-1)
#ax1.plot(ti,I0M,color=col[3],linestyle="dashed",linewidth=lw)
#ax11.plot(ti,I0P,color=col[3],label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$",linewidth=lw-1)
#ax11.plot(ti,I0M,color=col[3],linestyle="dashed",linewidth=lw-1)
ax1.set_yticks([0,0.2,0.4,0.6])
ax1.set_xlim(0,50)

#Is = I0P+I0M
#Id = I0P-I0M
#SP = 100 * Id ./ Is
#ax2.plot(ti,SP,color=col[3],linewidth=lw-1,label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$")
#ax22.plot(ti,SP,color=col[3],linewidth=lw-1,label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$")
#
#
#dat0P = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
#dat0M = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat")[2:end,:]
#
#l0P = length(dat0P[:,1])
#l0M = length(dat0M[:,1])
#
#lD = min(l0P,l0M)
#dat0P = dat0P[1:lD,:]
#dat0M = dat0M[1:lD,:]
#
#ti = dat0P[:,1]*fac_t
#I0P = (dat0P[:,2]-dat0P[:,3])*fac_c/2
#I0M = (dat0M[:,2]-dat0M[:,3])*fac_c/2
#
#ax1.plot(ti,I0P,color=col[4],label=L"$N=3\times 6, \ \lambda_1 = 0$",linewidth=lw-1)
#ax11.plot(ti,I0P,color=col[4],label=L"$N=3\times 6, \ \lambda_1 = 0$",linewidth=lw-1)
#ax1.plot(ti,I0M,color=col[4],linestyle="dashed",linewidth=lw)
#ax11.plot(ti,I0M,color=col[4],linestyle="dashed",linewidth=lw-1)
#ax1.set_yticks([0,0.2,0.4,0.6])
#
#Is = I0P+I0M
#Id = I0P-I0M
#SP = 100 * Id ./ Is
#ax2.plot(ti,SP,color=col[4],linewidth=lw-1,label=L"$N=3\times 6, \ \lambda_1 = 0$")
#ax22.plot(ti,SP,color=col[4],linewidth=lw-1,label=L"$N=3\times 6, \ \lambda_1 = 0$")
#
#ax1.set_ylabel("current "*L"$I\,[\mu\mathrm{A}]$",labelpad=15)	
#ax1.set_xticklabels([])
#ax1.set_xlim(0,50)
ax1.annotate("(a)",(0.02,0.86),xycoords="axes fraction",fontsize=fs)
#ax11.set_xlim(47,50)
#ax11.tick_params(axis="both", which="major", labelsize=fs-6)

ax2.set_ylabel("spin selectivity "*L"S\,[\%]")
ax2.set_xlim(0,50)
ax2.set_xticklabels([])
#ax2.set_yscale("log")
#ax2.set_ylim(1e-7,1e2)
ax2.annotate("(b)",(0.02,0.86),xycoords="axes fraction",fontsize=fs)
ax2.set_ylim(-80,100)

ax22.set_xlim(47,50)
ax22.set_yticks([-80,-40,0,40,80])
ax22.tick_params(axis="both", which="major", labelsize=fs-6)
sufP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat"
sufM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat"
for (p,c,l) in zip(pre,col,lab)
        datXP = readdlm("$(p)$(sufP)")[2:end,:]
        datXM = readdlm("$(p)$(sufM)")[2:end,:]
        lXP = length(datXP[:,1])
        lXM = length(datXM[:,1])
        lf = min(lXP,lXM)
	datXP = datXP[1:lf,:]
	datXM = datXM[1:lf,:]
        tim = datXP[:,1]*fac_t
	XP = datXP[:,2]+datXP[:,4]+datXP[:,6]+datXP[:,8]
	XM = datXM[:,2]+datXM[:,4]+datXM[:,6]+datXM[:,8]
	dX = XP-XM
        ax3.plot(tim,dX,color=c,linewidth=lw-1,label=l)
end

#dat0P = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
#dat0M = readdlm("L3x6/L18_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
#lXP = length(dat0P[:,1])
#lXM = length(dat0M[:,1])
#lf = min(lXP,lXM)
#dat0P = dat0P[1:lf,:]
#dat0M = dat0M[1:lf,:]
#ti = dat0P[:,1]*fac_t
#XP = dat0P[:,2]+dat0P[:,4]+dat0P[:,6]+dat0P[:,8]
#XM = dat0M[:,2]+dat0M[:,4]+dat0M[:,6]+dat0M[:,8]
#dX = XP-XM
#ax3.plot(ti,dX,color=col[3],linewidth=lw-1,label=L"$N=3\times 6, \ t_1=\lambda_1 = 0$")
#
#dat0P = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
#dat0M = readdlm("L3x6/L18_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")[2:end,:]
#lXP = length(dat0P[:,1])
#lXM = length(dat0M[:,1])
#lf = min(lXP,lXM)
#dat0P = dat0P[1:lf,:]
#dat0M = dat0M[1:lf,:]
#ti = dat0P[:,1]*fac_t
#XP = dat0P[:,2]+dat0P[:,4]+dat0P[:,6]+dat0P[:,8]
#XM = dat0M[:,2]+dat0M[:,4]+dat0M[:,6]+dat0M[:,8]
#dX = XP-XM
#ax3.plot(ti,dX,color=col[4],linewidth=lw-1,label=L"$N=3\times 6, \ \lambda_1 = 0$")
#

ax3.set_xlim(0,50)
ax3.set_ylim(-10,10)
ax3.set_ylabel("displacement polar. "*L"$S_X$")
ax3.set_xlabel("time "*L"$t\,[\mathrm{ps}]$")
ax3.annotate("(c)",(0.02,0.86),xycoords="axes fraction",fontsize=fs)

ax3.legend(loc="lower left",ncol=2)
plt.savefig("Cur_sp_ph_singleTrag_init_random_Nvar.pdf")
	
	
	
	






