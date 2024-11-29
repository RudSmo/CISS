using PyPlot
using DelimitedFiles
using FnMatch

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 13
lw = 4
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
rc("legend",fontsize=fs-5)
px = 0.16
py = 0.59
sx = 0.81
sx2= 0.35
sy = 0.38
dy = 0.1
dx = 0.08
fig = plt.figure(figsize=(4,6))

ax1=plt.axes([px,py,sx,sy])

sufP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"
sufM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"


pre = ["L1x5_init_random/L5"]
lab = [L"$N=1\times 5$",L"$N=3\times 6$"]
col = ["tab:blue","tab:orange","tab:green","tab:red"]
LW = [lw-1,lw+2,lw-1,lw-2]
for (p,c,l,lws) in zip(pre,col,lab,LW)
	datP = readdlm("$(p)$(sufP)")[2:end,:]
	datM = readdlm("$(p)$(sufM)")[2:end,:]
	lP = length(datP[:,1])
	lM = length(datM[:,1])
	lD = min(lP,lM)
	ti = datP[1:lD,1]*fac_t
	IP = (datP[1:lD,2]-datP[1:lD,3])*fac_c/2
	IM = (datM[1:lD,2]-datM[1:lD,3])*fac_c/2
	
	Is = IP+IM
	Id = IP-IM
	SP = 100 * Id ./ Is
	ax1.plot(ti,SP,color=c,label=l,linewidth=lws)
end


ax1.set_ylabel("spin selectivity "*L"S\,[\%]")
ax1.set_xlabel("time "*L"$t\, [\mathrm{ps}]$")
ax1.set_xlim(0,50)
ax1.set_ylim(-80,100)

plt.savefig("Cur_sp_ph_singleTrag_rand_init_L1x5.pdf")
	
	
	
	






