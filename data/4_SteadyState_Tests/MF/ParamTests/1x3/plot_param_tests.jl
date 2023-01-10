using LinearAlgebra
using PyPlot
using DelimitedFiles
using HDF5

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 5
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


fig = plt.figure(figsize = (10,8))
px = 0.075
px2 = 0.57
py2 = 0.7
py = 0.54
sx = 0.91
sx2 = 0.4
dy = 0.03
dy1 = dy-0.04
sy = 0.44
sy2 = 0.15

fac_c = 1e6 * 2.8494561
fac_t = 1/(4.135668*1000)
fac_s = 1e3

ax1 = plt.axes([px,py,sx,sy])
ax1.set_xticklabels([])

ax2 = plt.axes([px2,py2+0.02,sx2,sy2])
ax3 = plt.axes([px,py-sy-dy,sx,sy])
ax4 = plt.axes([px2,py2-sy-dy-0.09,sx2,sy2])


labels = [L"$t_1 = 4\,\mathrm{meV}$, $\lambda_0 = 0\,\mathrm{meV}$, $\lambda_1 = 0\,\mathrm{meV}$",
	  L"$t_1 = 4\,\mathrm{meV}$, $\lambda_0 = 1\,\mathrm{meV}$, $\lambda_1 = 0\,\mathrm{meV}$",
	  L"$t_1 = 4\,\mathrm{meV}$, $\lambda_0 = 1\,\mathrm{meV}$, $\lambda_1 = 0.1\,\mathrm{meV}$",
	  L"$t_1 = 0\,\mathrm{meV}$, $\lambda_0 = 1\,\mathrm{meV}$, $\lambda_1 = 0.1\,\mathrm{meV}$",
	  L"$t_1 = 0\,\mathrm{meV}$, $\lambda_0 = 1\,\mathrm{meV}$, $\lambda_1 = 0\,\mathrm{meV}$"]
colors = ["blue","green","red","turquoise","purple"]

ax = [[ax1,ax2],[ax3,ax4]]
fac_c = [fac_c,fac_c*1e-3]
for (V,ind) in zip([0.0,0.48],[1,2])
	folders = ["V$(V)/tPh_noSOC_noSOCPh/L3_tPh0.004_SOC0.0_SOCPh0.0_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat",
           	   "V$(V)/tPh_SOC_noSOCPh/L3_tPh0.004_SOC0.001_SOCPh0.0_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat",
           	   "V$(V)/tPh_SOC_SOCPh/L3_tPh0.004_SOC0.001_SOCPh0.0001_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat",
           	   "V$(V)/notPh_SOC_SOCPh/L3_tPh0.0_SOC0.001_SOCPh0.0001_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat",
            	   "V$(V)/notPh_SOC_noSOCPh/L3_tPh0.0_SOC0.001_SOCPh0.0_V$(V)_w00.002_P0.5_X0_-10_dt0.05_1_currents.dat"]	
	
	for (k,folder,label,col) in zip(Array(0:length(folders)),folders,labels,colors)
		dat = readdlm(folder)
		ax[ind][1].plot(fac_t * dat[2:end,1],fac_c[ind] * (dat[2:end,2]-dat[2:end,3]),linewidth=lw-0.7*k,color=col,label=label)
		ax[ind][2].plot(fac_t * dat[2:end,1],fac_c[ind] * (dat[2:end,2]-dat[2:end,3]),linewidth=lw-0.7*k,color=col,label=label)
	end
end

ax1.legend(loc="upper left",fontsize=fs-3)
ax3.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax1.set_ylabel(L"$I_L-I_R\,[\mathrm{nA}]$",labelpad=-2)

annx = 0.85
anny = 0.9
ax1.annotate(L"$V=0.0\, \mathrm{V}$",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax1.annotate("(a)",(0.03,0.44),xycoords="axes fraction",fontsize=fs)
ax3.set_ylabel(L"$I_L-I_R\,[\mu\mathrm{A}]$")
ax1.set_xlim(-0.05,3)
ax3.set_xlim(-0.05,3)
ax3.annotate(L"$V=0.48\, \mathrm{V}$",(annx,anny),xycoords="axes fraction",fontsize=fs)
ax3.annotate("(b)",(0.03,0.9),xycoords="axes fraction",fontsize=fs)
ax2.set_xlim(0.4,3)
ax4.set_xlim(0.4,3)
ax2.set_ylim(-0.01,0.01)
ax4.set_ylim(11,13.3)

		
plt.savefig("MF_1x3_ParamTests.pdf")


