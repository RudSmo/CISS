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

px = 0.12
py = 0.55
sx = 0.84
sy = 0.4
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,5))

Ls = [2,3]
t_E_0 = 1e3 * 0.2433946
t_E_1 = 1/(41.35668)
t_E_2 = 1/(0.002*41.35668)
I_E_1 = 1e3 * 0.2433946
I_E_2 = 0.002 * 1e3 * 0.2433946
t_T = 1e-3/2.418884326505E-2
ylim = [[1.1,1.3],[0.8,1.3],[0.3,0.75]]
#ylim = [[1.,1.4],[1.,1.4]]
for i in 1:length(Ls)
	k = i-1
	ax = plt.axes([px,py-k*(sy+dy),sx,sy])
#	TTHEOM = readdlm("L$(Ls[i])/TTHEOM/"*fnfilter(readdir("L$(Ls[i])/TTHEOM/"),"current*")[1])
	TAMEN= readdlm("L$(Ls[i])/TAMEN/res.txt")
#	EHRENFEST= readdlm("L$(Ls[i])/ehrenfest/currents.dat") 

	ax.plot(TAMEN[:,1],TAMEN[:,4],color="orange",label="maximal",linewidth=lw)
	ax.plot(TAMEN[:,1],TAMEN[:,3],color="green",label="mean",linewidth=lw)
#	ax.plot(TTHEOM[2:end,1]*t_E_1,(TTHEOM[2:end,2]-TTHEOM[2:end,3])*I_E_1,color="green",label="TTHEOM",linewidth=lw,linestyle="dashed")
	if i == 1
		t_E = t_E_1 
		I_E = I_E_1
	else	
		#t_E = t_E_1
		t_E = t_E_2
		I_E = I_E_2
	end
#	ax.plot(EHRENFEST[2:end,1],(EHRENFEST[2:end,2]-EHRENFEST[2:end,3])*I_E*2,color="turquoise",label="Ehrenfest",linewidth=lw,linestyle="dashed")	
	ax.annotate(L"$N=$"*"$(Ls[i])",(0.75,0.05),xycoords="axes fraction",fontsize=fs)
#	ax.set_xlim(0,40)
#	ax.set_ylim(ylim[i][1],ylim[i][2])
#	if i == 2
#		ax.legend()
#		ax.annotate(L"$V=0.48$"*" V", (0.05,0.05),xycoords="axes fraction",fontsize=fs)
#	end
	ax.set_ylabel("rank")
	if i == 1
		ax.set_xticklabels([])
		ax.legend(loc="center right")
	else
		ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
	end
end
plt.savefig("Peierls_LVar_V60_tamen_rank.pdf")



	


