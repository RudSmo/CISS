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
py = 0.13
sx = 0.84
sy = 0.8
dy = 0.05
dx = 0.01

lambda = 0.4
fig = plt.figure(figsize=(7,4))
Vs = [10]

t_E = 1/(41.35668)
I_E = 1e3 * 0.2433946
t_T = 1e-3/2.418884326505E-2
ylim = [[0.3,0.7]]
for i in 1:length(Vs)
	k = i-1
	ax = plt.axes([px,py-k*(sy+dy),sx,sy])
	TTHEOM = readdlm("V$(Vs[i])/TTHEOM/"*fnfilter(readdir("V$(Vs[i])/TTHEOM/"),"current*")[1])
	TAMEN00= readdlm("V$(Vs[i])/TAMEN/tol_e-2/curr00.txt")
	TAMEN01= readdlm("V$(Vs[i])/TAMEN/tol_e-2/curr01.txt")
	TAMEN_1  = zeros(length(TAMEN00[:,1]),2)
	TAMEN_1[:,1] = TAMEN00[:,1]
	TAMEN_1[:,2] = TAMEN00[:,2] - TAMEN01[:,2]

	ax.plot(TAMEN_1[:,1],TAMEN_1[:,2]*I_E,color="blue",label="TAMEN, tol="*L"$10^{-2}$",linewidth=lw)

	TAMEN00= readdlm("V$(Vs[i])/TAMEN/tol_e-3/curr00.txt")
        TAMEN01= readdlm("V$(Vs[i])/TAMEN/tol_e-3/curr01.txt")
        TAMEN_2  = zeros(length(TAMEN00[:,1]),2)
        TAMEN_2[:,1] = TAMEN00[:,1]
        TAMEN_2[:,2] = TAMEN00[:,2] - TAMEN01[:,2]	
	ax.plot(TAMEN_2[:,1],TAMEN_2[:,2]*I_E,color="orange",label="TAMEN, tol="*L"$10^{-3}$",linewidth=lw)

	ax.plot(TTHEOM[2:end,1]*t_E,(TTHEOM[2:end,2]-TTHEOM[2:end,3])*I_E,color="red",label="TTHEOM",linewidth=lw,linestyle="dashed")
	ax.annotate(L"$V=$"*"$(Vs[i]*0.1) V",(0.75,0.05),xycoords="axes fraction",fontsize=fs)
	ax.set_xlim(0,50)
	ax.set_ylim(ylim[i][1],ylim[i][2])
	#if i == 1
		ax.legend()
		ax.annotate(L"$\lambda=$"*"$(lambda) eV", (0.05,0.05),xycoords="axes fraction",fontsize=fs)
	#end
	ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
	#if i != 3
	#	ax.set_xticklabels([])
	#else
		ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
	#end
end
plt.savefig("Holstein_lambda_$(lambda)_V10_conv.pdf")



	


