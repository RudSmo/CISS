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

fig = plt.figure(figsize=(10,7))

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

ax1.set_xticklabels([])

dat = readdlm("L8_tPh0.0_SOC0.001_SOCPh0.0_V0.0_w00.002_P0.5_X0_-1.0_dt0.05_1_electrons.dat")
colors = plt.cm.get_cmap("rainbow",8)
colors = [colors(i) for i in 0:1/7:1]
for (i,c) in zip(1:8,colors)
	ax1.plot(dat[2:end,1]*fac_t,dat[2:end,2*i],color=c,linewidth=lw,label="Site $(i)")
	ax2.plot(dat[2:end,1]*fac_t,dat[2:end,2*i+1],color=c,linewidth=lw)
end
ax1.legend(loc="best",ncol=4)
ax1.set_ylim(0.94,1.01)
ax1.set_yticks(0.94:0.01:1.0)
ax1.set_yticklabels([string(i) for i in 0.94:0.01:1.0])
ax1.set_xlim(0,1.8)
ax2.set_xlim(0,1.8)
ax2.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax1.set_ylabel(L"$n_{\uparrow}$")
ax2.set_ylabel(L"$n_{\downarrow}$")
plt.savefig("MF_noPh_2x4_population_dynamics_V0.pdf")

