using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
lw = 2
fac_c = 0.002* 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.09
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(9,4.5))
ax = plt.axes([px,py,sx,sy])

fac_c = 0.002* 1e3 * 0.2433946
colors = ["k","red","blue","green","orange","purple","magenta","turquoise"]
syssize = [3,4,6,8,9,10,12,12]
sizes = ["L1x3","L2x2","L2x3","L2x4","L3x3","L2x5","L2x6","L3x4"]
labels = [L"$L=1\times3$",L"$L=2\times 2$",L"$L=2\times 3$",L"$L=2\times 4$",L"$L=3\times 3$",L"$L=2\times 5$",L"$L=2\times 6$",L"$L=3\times 4$"]
Ns = [71,10,72,11,11,11,11,11]
ms = 8
fac_V = 0.2 * 0.04
for (sys,size,lab,col,ns) in zip(syssize,sizes,labels,colors,Ns)
	datP = readdlm("$(size)/Phonons_steadyState_$(ns)_P.dat")
	datM = readdlm("$(size)/Phonons_steadyState_$(ns)_M.dat")	
	ax.plot(datP[:,1]*fac_V,datP[:,2],color=col,label=lab,linewidth=lw)
	ax.plot(datM[:,1]*fac_V,datM[:,2],color=col,linestyle="dashed",linewidth=lw)
end



ax.set_xlabel(L"$V\,[\mathrm{V}]$")
ax.set_ylabel(L"$X_1$")
ax.legend(loc="lower left")
#plt.savefig("PhX1_V_sizeVar_P$(P)_VVar.png",dpi=200)
plt.savefig("PhX1_V_sizeVar_PVar_VVar.pdf")
