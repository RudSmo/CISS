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

rc("axes",labelsize=fs-2)
rc("xtick",labelsize=fs-4)
rc("ytick",labelsize=fs-4)

px = 0.18
py = 0.59
sx = 0.79
sx2= 0.35
sy = 0.38
dy = 0.1
dx = 0.08
fig = plt.figure(figsize=(4,6))

ax1=plt.axes([px,py,sx,sy])
ax11=plt.axes([px,py-sy-dy,sx,sy])
#ax12=plt.axes([px+sx2+dx,py-sy-dy,sx2,sy])
datP=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")
datM=readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")

start = 2#10000
ende = length(datP[:,1])
XtP = (datP[start:ende,2] .+ datP[start:ende,4] .+ datP[start:ende,6] .+ datP[start:ende,8]) ./4
PtP = (datP[start:ende,3] .+ datP[start:ende,5] .+ datP[start:ende,7] .+ datP[start:ende,9]) ./4   
XtM = (datM[start:ende,2] .+ datM[start:ende,4] .+ datM[start:ende,6] .+ datM[start:ende,8]) ./4
PtM = (datM[start:ende,3] .+ datM[start:ende,5] .+ datM[start:ende,7] .+ datM[start:ende,9]) ./4
println(datP[ende,1]*fac_t)
ax1.plot(XtP,PtP,color="red",label=L"$s=\uparrow$")
ax1.plot(XtM,PtM,color="blue",linestyle="dashed",label=L"$s=\downarrow$")
ax1.set_xlabel(L"$\bar{X}$")
ax1.set_ylabel(L"$\bar{P}$",labelpad=-5)
ax1.legend(loc="best",fontsize=fs-5,frameon=false)

dX = XtP-XtM
dP = PtP-PtM

ax11.plot(dX,dP,color="k")
ax11.set_xlabel(L"$\bar{X}_+-\bar{X}_-$")
ax11.set_ylabel(L"$\bar{P}_+-\bar{P}_-$",labelpad=-2)
ax1.set_xlim(-50,60)
ax1.set_ylim(-50,66)
#ax11.set_xlim(-0.5,0.5)
#ax11.set_ylim(-0.5,0.5)

#plt.show()

plt.savefig("PhaseSpace_L5_init_X0P0_full.pdf")
	
	
	
	






