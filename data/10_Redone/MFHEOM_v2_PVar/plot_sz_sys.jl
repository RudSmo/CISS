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
syssize = [1:3,1:4,1:6,1:8,1:9,1:10,1:12,1:12]
sizes = ["L1x3","L2x2","L2x3","L2x4","L3x3","L2x5","L2x6","L3x4"]
labels = [L"$L=1\times3$",L"$L=2\times 2$",L"$L=2\times 3$",L"$L=2\times 4$",L"$L=3\times 3$",L"$L=2\times 5$",L"$L=2\times 6$",L"$L=3\times 4$"]
Ns = [71,10,72,11,11,11,11,11]
ms = 8
fac_V = 0.2 * 0.04
V = 65
for (sys,size,lab,col,ns) in zip(syssize,sizes,labels,colors,Ns)
	datP = readdlm("$(size)/Electrons_steadyState_$(ns)_P.dat")
	datM = readdlm("$(size)/Electrons_steadyState_$(ns)_M.dat")	
	indP,indM = 1,1
	for i in 1:length(datP[:,1])	
		if abs(datP[i,1] - V) < 1e-5
			indP = i
		end
	end
	for i in 1:length(datM[:,1])
                if abs(datM[i,1] - V) < 1e-5
                        indM = i
                end
        end	
	ax.plot(sys,[datP[indP,j]-datP[indP,j+1] for j in 2:2:length(datP[1,:])],color=col,label=lab,linewidth=lw)
	ax.plot(sys,[datM[indM,j]-datM[indM,j+1] for j in 2:2:length(datM[1,:])],color=col,linestyle="dashed",linewidth=lw)
end



ax.set_xlabel("site")
ax.set_ylabel(L"$s_i^z$")
ax.legend(loc="lower right",ncol=2)
ax.set_xlim(0.5,12.5)
#plt.savefig("Sz1_V_sizeVar_P$(P)_VVar.png",dpi=200)
plt.savefig("Sz_space_V$(V)_sizeVar_PVar.pdf")
