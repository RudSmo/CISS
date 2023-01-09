using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
lw = 3
fac_t = 1/(4.135668*1000) 


px = 0.08
py = 0.54
sx = 0.89
sy = 0.44
dy = 0.01

px2 = 0.5
py2 = 0.63
sx2 = 0.4
sy2 = 0.25

fig = plt.figure(figsize=(10,7))

ax1 = plt.axes([px,py,sx,sy])
ax11 = plt.axes([px2,py2,sx2,sy2])
ax2 = plt.axes([px,py-sy-dy,sx,sy])
ax22 = plt.axes([px2,py2-sy-dy,sx2,sy2])

ax1.set_xticklabels([])

dat1 = readdlm("L8_tPh0.0_SOC0.001_SOCPh0.0_V0.0_w00.002_P0.0_X0_-1.0_dt0.05_1_electrons.dat")
dat2 = readdlm("L8_tPh0.0_SOC0.0_SOCPh0.0_V0.0_w00.002_P0.0_X0_-1.0_dt0.05_1_electrons.dat")

ax1.plot(dat1[2:end,1]*fac_t,dat1[2:end,2],linewidth=lw,color="k",label=L"$\lambda_0=0.04\,\mathrm{eV}$")
ax1.plot(dat2[2:end,1]*fac_t,dat2[2:end,2],linewidth=lw-1,color="red",label=L"$\lambda_0=0.0\,\mathrm{eV}$")
ax11.plot(dat1[2:end,1]*fac_t,dat1[2:end,2],linewidth=lw,color="k",label=L"$\lambda=0.04\,\mathrm{eV}$")
ax11.plot(dat2[2:end,1]*fac_t,dat2[2:end,2],linewidth=lw-1,color="red",label=L"$\lambda=0.0\,\mathrm{eV}$")
ax2.plot(dat1[2:end,1]*fac_t,dat1[2:end,3],linewidth=lw,color="k",label=L"$\lambda=0.04\,\mathrm{eV}$")
ax2.plot(dat2[2:end,1]*fac_t,dat2[2:end,3],linewidth=lw-1,color="red",label=L"$\lambda=0.0\,\mathrm{eV}$")
ax22.plot(dat1[2:end,1]*fac_t,dat1[2:end,3],linewidth=lw,color="k",label=L"$\lambda=0.04\,\mathrm{eV}$")
ax22.plot(dat2[2:end,1]*fac_t,dat2[2:end,3],linewidth=lw-1,color="red",label=L"$\lambda=0.0\,\mathrm{eV}$")


ax1.legend(loc="best",ncol=4)
ax1.set_ylim(0.96,0.995)
ax1.set_yticks(0.94:0.01:1.0)
ax1.set_yticklabels([string(i) for i in 0.94:0.01:1.0])
ax1.set_xlim(0,1.8)
ax1.set_ylabel(L"$n_{1\uparrow}$")

ax11.set_xlim(0.5,1.8)
ax11.set_ylim(0.992,0.9935)
ax11.set_yticks(0.992:0.0005:0.9935)
ax11.set_yticklabels([string(i) for i in 0.9920:0.0005:0.9935])

ax2.set_xlim(0,1.8)
ax2.set_yticks(0.96:0.01:0.99)
ax2.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax2.set_ylabel(L"$n_{1\downarrow}$")
ax22.set_xlim(0.5,1.8)
ax22.set_ylim(0.992,0.9935)
ax22.set_yticks(0.992:0.0005:0.9935)
ax22.set_yticklabels([string(i) for i in 0.9920:0.0005:0.9935])

plt.savefig("MF_noPh_2x4_population_dynamics_V0_P0_LambdaVar.pdf")

