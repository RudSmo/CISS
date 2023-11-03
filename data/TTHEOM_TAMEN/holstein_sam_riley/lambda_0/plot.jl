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
py = 0.14
sx = 0.84
sy = 0.8
dy = 0.05
dx = 0.01

lambda = 0.0
fig = plt.figure(figsize=(7,4))
Vs = [1]

t_E = 1/(41.35668)
I_E = 1e3 * 0.2433946
t_T = 1/2.418884326505E-2
#ylim = [[-0.01,0.01],[0.6,1.3],[1.5,2.1]]
for i in 1:length(Vs)
	k = i-1
	ax = plt.axes([px,py-k*(sy+dy),sx,sy])
	TTHEOM = readdlm("V$(Vs[i])/TTHEOM/"*fnfilter(readdir("V$(Vs[i])/TTHEOM/"),"current*")[1])
	TAMEN00= readdlm("V$(Vs[i])/TAMEN/curr00.txt")
	TAMEN01= readdlm("V$(Vs[i])/TAMEN/curr01.txt")
	TAMEN  = zeros(length(TAMEN00[:,1]),2)
	TAMEN[:,1] = TAMEN00[:,1]
	TAMEN[:,2] = TAMEN00[:,2] - TAMEN01[:,2]
	EHREN  = readdlm("V$(Vs[i])/Ehrenfest/data0/L1_gPhU0.0_gPhD0.0_V1.0_w00.2_dt0.02_currents_avg.dat")

	ax.plot(TAMEN[:,1],TAMEN[:,2]*I_E,color="blue",label="TAMEN",linewidth=lw)
	ax.plot(TTHEOM[2:end,1]*t_E,(TTHEOM[2:end,2]-TTHEOM[2:end,3])*I_E,color="red",label="TTHEOM",linewidth=lw)
	ax.plot(EHREN[2:end,1]*t_E,(EHREN[2:end,2]-EHREN[2:end,3])*I_E*2,color="k",label="Ehrenfest",linewidth=lw)
	ax.annotate(L"$V=$"*"$(Vs[i]) V",(0.75,0.05),xycoords="axes fraction",fontsize=fs)
	ax.set_xlim(0,5)
	ax.set_ylim(2.2,2.8)
#	ax.set_ylim(ylim[i][1],ylim[i][2])
	if i == 1
		ax.legend()
		ax.annotate(L"$\lambda=$"*"$(lambda) eV", (0.2,0.05),xycoords="axes fraction",fontsize=fs)
	end
	ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
	#if i != 3
#		ax.set_xticklabels([])
	#else
	ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
	#end
end
plt.savefig("Holstein_lambda_$(lambda)_VVar.pdf")



	


