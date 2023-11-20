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

px = 0.14
py = 0.72
sx = 0.84
sy = 0.26
dy = 0.05
dx = 0.01

lambda = 0.3
fig = plt.figure(figsize=(7,6))
Vs = [5,10]

t_E = 1/(41.35668)
#I_E = 1.0
I_E = 1e3 * 0.2433946
t_T = 1e-3/2.418884326505E-2
ylim = [[0.0,1.1],[0.0,2.3]]
for i in 1:length(Vs)
	k = i-1
	ax = plt.axes([px,py-k*(sy+dy),sx,sy])
	TTHEOM = readdlm("V$(Vs[i])/TTHEOM/"*fnfilter(readdir("V$(Vs[i])/TTHEOM/"),"n_vib*")[1])
	if i == 1
		TTHEOM_inds = sortperm(TTHEOM[2:end,1])
                TTHEOM = TTHEOM[2:end,:][TTHEOM_inds,:][1:10000,:]
	end
#	TAMEN= readdlm("V$(Vs[i])/TAMEN/tol_1e-3/curr.txt")
	SAM=readdlm("../propagation_check/nvib_lambda_300_voltage_$(0.1*Vs[i])V.txt")
#	EHREN  = readdlm("V$(Vs[i])/Ehrenfest/currents_avg_n33000.dat")

#	ax.plot(TAMEN[:,1],(TAMEN[:,2]-TAMEN[:,3])*I_E,color="orange",label="TAMEN",linewidth=lw)
	ax.plot(TTHEOM[2:end,1]*t_E,TTHEOM[2:end,2],color="green",label="TTHEOM",linewidth=lw)
	ax.plot(SAM[:,1]*t_E,SAM[:,2], color="orange", label="Sams HEOM", linewidth=lw, linestyle="dashed")
#	ax.plot(EHREN[:,1]*t_E,(EHREN[:,2]-EHREN[:,3])*I_E*2,color="k",label="Ehrenfest",linewidth=lw)
	ax.annotate(L"$V=$"*"$(Vs[i]*0.1) V",(0.75,0.05),xycoords="axes fraction",fontsize=fs)
	ax.set_xlim(0,50)
	ax.set_ylim(ylim[i][1],ylim[i][2])
	if i == 1
		ax.legend()
		ax.annotate(L"$\lambda=$"*"$(lambda) eV", (0.05,0.05),xycoords="axes fraction",fontsize=fs)
	end
	ax.set_ylabel(L"$\langle a^\dagger a\rangle$")
	if i != 2
		ax.set_xticklabels([])
	else
		ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
	end
end
plt.savefig("Holstein_lambda_$(lambda)_Nvib_VVar.pdf")



	


