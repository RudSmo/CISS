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
px = 0.08
py = 0.72
sx1 = 0.285
dx1 = 0.02
sx = 0.43
dx = 0.038
dy = 0.05
sy = 0.275
ax11 = plt.axes([px,py,sx1,sy])
ax12 = plt.axes([px+sx1+dx1,py,sx1,sy]) 
ax13 = plt.axes([px+2*(sx1+dx1),py,sx1,sy])
ax2 = plt.axes([px,py-sy-dy,2*sx+dx,sy])
ax31 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax32 = plt.axes([px+sx+dx,py-2*(sy+dy),sx,sy])

folds = ["X0_0.0","X0_-3.0","X0_-10.5","X0_-25.0"]
fold_str = [L"$0$",L"$-3$",L"$-10.5$",L"$-25$"]
a = 10000
lw = 2
fac = 1e0
fac_c = 1e0
for (f,fstr) in zip(folds,fold_str)
	dat_pop = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V12.0_$(f)_dt0.002_electrons.dat")
	ax11.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,2]-dat_pop[2:end-a,3])*fac,linewidth=lw)
	ax12.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,4]-dat_pop[2:end-a,5])*fac,linewidth=lw)
	ax13.plot(dat_pop[2:end-a,1],(dat_pop[2:end-a,6]-dat_pop[2:end-a,7])*fac,linewidth=lw)
	dat_cur = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V12.0_$(f)_dt0.002_currents.dat")
        ax2.plot(dat_cur[2:end,1],(dat_cur[2:end,2]-dat_cur[2:end,3])*fac_c,linewidth=lw,label=L"$X_0=$"*fstr)

	dat_ph = readdlm("$(f)/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V12.0_$(f)_dt0.002_phonons.dat")
	ax31.plot(dat_ph[2:end,2],dat_ph[2:end,3],linewidth=lw)
	ax32.plot(dat_ph[2:end,4],dat_ph[2:end,5],linewidth=lw)
end

datSS_E = readdlm("../../../../../noPh/MF-HEOM/1x3/P0.5/Electrons_steadyState.dat")[121,2:7]
datSS_C = readdlm("../../../../../noPh/MF-HEOM/1x3/P0.5/Currents_steadyState.dat")[121,2:3]

lw = 2

ax11.plot(Array(0:4900),fac*(datSS_E[1]-datSS_E[2])*ones(length(Array(0:4900))),color="orange",linewidth=2,linestyle="dashed")
ax11.plot([4900],fac*(datSS_E[1]-datSS_E[2]),color="orange",marker="o",markersize=8,linewidth=0)
ax12.plot(Array(0:4900),fac*(datSS_E[3]-datSS_E[4])*ones(length(Array(0:4900))),color="orange",linewidth=2,linestyle="dashed")
ax12.plot([4900],fac*(datSS_E[3]-datSS_E[4]),color="orange",marker="o",markersize=8,linewidth=0)
ax13.plot(Array(0:4900),fac*(datSS_E[5]-datSS_E[6])*ones(length(Array(0:4900))),color="orange",linewidth=2,linestyle="dashed")
ax13.plot([4900],fac*(datSS_E[5]-datSS_E[6]),color="orange",marker="o",markersize=8,linewidth=0)

ax11.set_ylabel(L"$s^z$")
ax11.set_ylim(-0.03,0.17)
ax11.set_xticks([0,1000,2000,3000,4000])
ax11.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax11.annotate("(a)",(0.05,0.9),xycoords="axes fraction",fontsize=fs)
ax11.annotate("site 1",(0.76,0.89),xycoords="axes fraction",fontsize=fs)

ax12.set_ylim(-0.03,0.17)
ax12.set_yticklabels([])
ax12.set_xticks([0,1000,2000,3000,4000])
ax12.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax12.annotate("(b)",(0.05,0.9),xycoords="axes fraction",fontsize=fs)
ax12.annotate("site 2",(0.76,0.89),xycoords="axes fraction",fontsize=fs)

ax13.set_ylim(-0.03,0.17)
ax13.set_yticklabels([])
ax13.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=13,labelpad=-2)
ax13.set_xlim(0,5000)
ax13.annotate("(c)",(0.05,0.9),xycoords="axes fraction",fontsize=fs)
ax13.annotate("site 3",(0.76,0.89),xycoords="axes fraction",fontsize=fs)

ax2.plot(Array(0:7000),fac_c*(datSS_C[1]-datSS_C[2])*ones(length(Array(0:7000))),color="orange",linewidth=2,linestyle="dashed",label=L"Steady state $t_1 = \lambda_1 = 0$")
ax2.plot([6950],fac_c*(datSS_C[1]-datSS_C[2]),color="orange",marker="o",markersize=8,linewidth=0)
ax2.plot([-10000],[-1],marker="x",color="red",linewidth=0,markersize="8",label=L"$\xi^\mathrm{adiab.}_\mathrm{ss}$ for $X_0 = -10.5$")
ax2.set_ylim(-0.01,0.2)
ax2.set_xlim(0,7000)
ax2.set_ylabel(L"$I_L-I_R \, [t_0]$")
ax2.set_xlabel(L"$t\, [t_0^{-1}]$",fontsize=15,labelpad=-2)
ax2.annotate("(d)",(0.02,0.85),xycoords="axes fraction",fontsize=fs)
ax2.legend(loc="upper right",ncol=2,numpoints=1)


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.]
    return a[1]*sx + a[2]*sy + a[3]*sz
end

function absv(a)
    return sqrt(dot(a,a))
end

function UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1)
    cj = c*(j-1)/(N-1)
    return [a*cos(phij),a*sin(phij),cj]
end

function VSigma(tso,a,c,k,num_laps,N)
        rj = UnitVector(a,c,k,num_laps,N)
        rj1 = UnitVector(a,c,k+1,num_laps,N)
        rj2 = UnitVector(a,c,k+2,num_laps,N)
        rjrj1 = rj-rj1
        rjrj2 = rj-rj2

        dj1 = rjrj1/absv(rjrj1)
        dj2 = rjrj2/absv(rjrj2)

        vj = cross(dj1,dj2)

        Vj = 1im*tso*VDot(vj)
        return Vj
end

tso = 1/40
tph = 0.1
w0 = 0.05
rho = h5open("../../../../../../../src/MeanField_v2_weighting_noPhDyn/NoPh_Coherences/L3_tPh0.1_SOC0.025_SOCPh0.0025000000000000005_V12.0_w00.05_P0.5_X0_-10.5_dt0.002_1_final_Quantities")["Rho"][]
t_rho_1 = -tph * sqrt(2*w0) * tr(rho[1:2,3:4])
vs1 = VSigma(tso,1.0,1.0,1,1,3)
l_rho_1 = tr(vs1 * rho[1:2,5:6])
X1 = 1/w0 * (t_rho_1 + l_rho_1)
X1 += conj(X1)
println(X1)


t_rho_2 = -tph * sqrt(2*w0) * tr(rho[3:4,5:6])
X2 = 1/w0 * t_rho_2
X2 += conj(X2)
println(X2)

ax31.annotate("(e)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax31.plot([X1],[0],color="red",marker="x",markersize=9,linewidth=3)
ax31.set_ylabel(L"$P_1$")
ax31.set_xlabel(L"$X_1$")
ax31.set_ylim(-12,12)
ax31.set_xlim(-40,30)

ax32.annotate("(f)",(0.06,0.85),xycoords="axes fraction",fontsize=fs)
ax32.plot([X2],[0],color="red",marker="x",markersize=9,linewidth=3)
ax32.set_yticklabels([])
ax32.set_ylim(-12,12)
ax32.set_xlim(-40,30)
ax32.set_xlabel(L"$X_2$")
ax32.set_ylabel(L"$P_2$",labelpad=-2)

plt.savefig("Pops_currents_phonons_1x3_V12.pdf")





