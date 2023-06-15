using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
num_sites = 3
num_laps = 2
N = num_laps * num_sites
V = 0.04 * 0.2 * 48
P = 0.5
Vs = ["V50","V60","V75"]
Ns = [72,71,72]
labels = [L"$V=\,0.4\mathrm{V}$",L"$V=\,0.48\mathrm{V}$",L"$V=\,0.6\mathrm{V}$"]
colors = ["blue","red","black"]
colors2 = ["turquoise","magenta","orange"]
ta = 25000
for (v,lab,col,col2,nS) in zip(Vs,labels,colors,colors2,Ns)
        dat = readdlm("$(v)/currents_avg_n$(nS)_P.dat")
	dat_s = readdlm("$(v)/currents_std_n$(nS)_P.dat")
	cur = (dat[:,2]-dat[:,3])*fac_c
	cur_s = sqrt.((dat_s[:,2]*fac_c).^2 + (dat_s[:,3]*fac_c).^2) / sqrt(nS)
        ax.plot(dat[:,1]*fac_t*1e-3,cur,linewidth=lw,label=lab,color=col)
	ax.fill_between(dat[:,1]*fac_t*1e-3,cur,y2=cur+cur_s,color=col,alpha=0.4)
	ax.fill_between(dat[:,1]*fac_t*1e-3,cur,y2=cur-cur_s,color=col,alpha=0.4)
	ax.plot(dat[ta:end,1]*fac_t*1e-3,cur[ta:end],linewidth=lw,color=col2)	
end
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$\langle I(t)\rangle\,[\mu\mathrm{A}]$")
ax.set_ylim(0.0,2.2)
ax.legend()


plt.savefig("L$(num_laps)x$(num_sites)_P$(P)_VVar_current.pdf")
