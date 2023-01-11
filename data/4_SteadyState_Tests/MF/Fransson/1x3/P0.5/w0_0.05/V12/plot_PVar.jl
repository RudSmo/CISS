using LinearAlgebra
using PyPlot
using DelimitedFiles
using HDF5

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

fig = plt.figure(figsize = (10,10))
px = 0.072
py = 0.71
sx1 = 0.283
dx1 = 0.03
sx = 0.435
dx = 0.038
dy = 0.05
sy = 0.275
ax11 = plt.axes([px,py,sx1,sy])
ax12 = plt.axes([px+sx1+dx1,py,sx1,sy]) 
ax13 = plt.axes([px+2*(sx1+dx1),py,sx1,sy])
ax2 = plt.axes([px,py-sy-dy,2*sx+dx,sy])
ax31 = plt.axes([px,py-2*(sy+dy),sx-0.015,sy])
ax32 = plt.axes([px+sx+dx+0.005,py-2*(sy+dy),sx-0.015,sy])

a = 0 
lw = 2
fac = 1
fac_c = 1e3 * 2.8494561
fac_t = 1/(4.135668*1000)

prefixes = ["","../../../P-0.5/w0_0.05/V12/"]
labels = [L"$P=0.5$",L"$P=-0.5$"]

for (pre,lab) in zip(prefixes,labels)
	dat_pop = readdlm("$(pre)electrons_avg_n10.dat")
	ax11.plot(dat_pop[2:end-a,1]*fac_t,(dat_pop[2:end-a,2]-dat_pop[2:end-a,3])*fac,linewidth=lw,label=lab)
	ax12.plot(dat_pop[2:end-a,1]*fac_t,(dat_pop[2:end-a,4]-dat_pop[2:end-a,5])*fac,linewidth=lw)
	ax13.plot(dat_pop[2:end-a,1]*fac_t,(dat_pop[2:end-a,6]-dat_pop[2:end-a,7])*fac,linewidth=lw)

	dat_ph = readdlm("$(pre)phonons_avg_n10.dat")
	ax31.plot(dat_ph[2:end,2],dat_ph[2:end,3],linewidth=lw)
	ax32.plot(dat_ph[2:end,4],dat_ph[2:end,5],linewidth=lw)
end
datSS_E = readdlm("../../../../../../noPh/V0.48/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.025_1_electrons.dat")[end,2:end]
datSS_C = readdlm("../../../../../../noPh/V0.48/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.025_1_currents.dat")[end,2:3]

dat_c1 = readdlm("$(prefixes[1])currents_avg_n10.dat")
cur1 = (dat_c1[2:end,2]-dat_c1[2:end,3])
dat_c2 = readdlm("$(prefixes[2])currents_avg_n10.dat")
cur2 = (dat_c2[2:end,2]-dat_c2[2:end,3])

cur_dif = cur1-cur2
cur_tot = cur1+cur2
mag_res = cur_dif ./ cur_tot * 100


ax2.plot(dat_c1[2:end,1]*fac_t, mag_res, linewidth=lw)

lw = 2

t_end = 60
#ax11.plot(Array(0:t_end),fac*(datSS_E[1]-datSS_E[2])*ones(length(Array(0:t_end))),color="orange",linewidth=2,linestyle="dashed",label=L"Steady state $t_1 = \lambda_1 = 0$")
#ax11.plot([t_end],fac*(datSS_E[1]-datSS_E[2]),color="orange",marker="o",markersize=8,linewidth=0)
#ax12.plot(Array(0:t_end),fac*(datSS_E[3]-datSS_E[4])*ones(length(Array(0:t_end))),color="orange",linewidth=2,linestyle="dashed")
#ax12.plot([t_end],fac*(datSS_E[3]-datSS_E[4]),color="orange",marker="o",markersize=8,linewidth=0)
#ax13.plot(Array(0:t_end),fac*(datSS_E[5]-datSS_E[6])*ones(length(Array(0:t_end))),color="orange",linewidth=2,linestyle="dashed")
#ax13.plot([t_end],fac*(datSS_E[5]-datSS_E[6]),color="orange",marker="o",markersize=8,linewidth=0)

ax11.set_ylabel(L"$s^z$")
ax11.set_ylim(-0.14,0.14)
ax11.set_xlim(0,t_end)
#ax11.set_xticks([0,1000,2000,3000,4000])
ax11.set_xlabel(L"$t\, [\mathrm{ps}]$",fontsize=13)#,labelpad=-1)
ax11.annotate("(a)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax11.annotate("site 1",(0.75,0.88),xycoords="axes fraction",fontsize=fs)
ax11.legend(loc="center right",numpoints=1,fontsize=12)

ax12.set_ylim(-0.14,0.14)
ax12.set_yticklabels([])
ax12.set_xlim(0,t_end)
#ax12.set_xticks([0,1000,2000,3000,4000])
ax12.set_xlabel(L"$t\, [\mathrm{ps}]$",fontsize=13)
ax12.annotate("(b)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax12.annotate("site 2",(0.75,0.88),xycoords="axes fraction",fontsize=fs)

ax13.set_ylim(-0.14,0.14)
ax13.set_xlim(0,t_end)
ax13.set_yticklabels([])
ax13.set_xlabel(L"$t\, [\mathrm{ps}]$",fontsize=13)#,labelpad=-2)
ax13.annotate("(c)",(0.1,0.85),xycoords="axes fraction",fontsize=fs)
ax13.annotate("site 3",(0.75,0.87),xycoords="axes fraction",fontsize=fs)

#ax2.plot(Array(0:t_end)*fac_t,fac_c*(datSS_C[1]-datSS_C[2])*ones(length(Array(0:t_end))),color="orange",linewidth=2,linestyle="dashed")
#ax2.plot([t_end],fac_c*(datSS_C[1]-datSS_C[2]),color="orange",marker="o",markersize=8,linewidth=0)
#ax2.set_ylim(-0.2,0.2)
ax2.set_xlim(0,t_end)
ax2.set_ylabel(L"$M_R = \frac{I_\uparrow-I_\downarrow}{I_\uparrow+I_\downarrow} \, [\%]$",labelpad=-2)
ax2.set_xlabel(L"$t\, [\mathrm{ps}]$",fontsize=15)#,labelpad=-2)
ax2.annotate("(d)",(0.02,0.85),xycoords="axes fraction",fontsize=fs)

ax31.annotate("(e)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax31.set_ylabel(L"$P_1$",labelpad=-1)
ax31.set_xlabel(L"$X_1$")

ax32.annotate("(f)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax32.set_yticklabels([])
ax32.set_xlabel(L"$X_2$")
ax32.set_ylabel(L"$P_2$",labelpad=-2)

plt.savefig("Pops_currents_phonons_1x3_Avg_V0.48_PComp.pdf")





