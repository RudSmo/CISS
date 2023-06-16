using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002 * 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002 * 4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,6))
ax = plt.axes([px,py,sx,sy])

#Ns =[[11,21,21,72,10,71],
#    [11,19,21,71,10,61],
#    [11,20,11,72,10,61]]
Ns = [[11,19,21,71,10,61]]

colors = ["k","red","blue","green","orange","purple"]
sizes = ["L1x3","L2x2","L2x3","L2x4","L3x3","L2x5"]
labels = [L"$L=1\times3$",L"$L=2\times 2$",L"$L=2\times 3$",L"$L=2\times 4$",L"$L=3\times 3$",L"$L=2\times 5$"]
V = [60]#[50,60,75]
P = 0.5
for (ns,v) in zip(Ns,V)
	for (nS,size,col,label) in zip(ns[end:-1:1],sizes,colors,labels)
		dat = readdlm("$(size)/V$(v)/currents_avg_n$(nS)_P.dat")
		ax.plot(dat[2:end,1]*fac_t*1e-3,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw,label=label,color=col)
	end
end
ax.annotate(L"$V=0.48\,\mathrm{V}$",(0.1,0.9),xycoords="axes fraction",fontsize=fs)		
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
ax.legend()

plt.savefig("Current_dyn_P$(P)_sysVar_V60.pdf")
