using PyPlot
using DelimitedFiles
using FnMatch

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.84
sx = 0.87
sy = 0.13
dy = 0.02

fig = plt.figure(figsize=(10,10*sqrt(2)))

ws = 0.002

fac_t_T = 1/(4.135668)
fac_t_E = 1/(ws*4.135668)
fac_c_T = 1e3 * 0.2433946/2
fac_c_E = 0.002* 1e3 * 0.2433946

annx = 0.85
anny = 0.05

Vs = [22,50,60,75,85,105]
yl = [[-0.05,0.05],[0.25,0.4],[0.45,0.7],[0.8,1.12],[0.9,1.3],[1.0,1.35]]
for (i,V) in enumerate(Vs)
	ax = plt.axes([px,py-(i-1)*(sy+dy),sx,sy])
	datE = readdlm("Ehrenfest/V$(V)/$(fnfilter(readdir("Ehrenfest/V$(V)"),"*currents*")[1])")
	datT = readdlm("TTHEOM/V$(V)/$(fnfilter(readdir("TTHEOM/V$(V)"),"*current*")[1])") 	
	datN = readdlm("noPh/V$(V)/$(fnfilter(readdir("noPh/V$(V)"),"*currents*")[1])")
	
	ax.plot(datE[2:end,1]*fac_t_E, (datE[2:end,2]-datE[2:end,3])*fac_c_E, color="blue",linewidth=lw,label="Ehrenfest")
	ax.plot(datT[2:end,1]*fac_t_T, (datT[2:end,2]-datT[2:end,3])*fac_c_T, color="green",linewidth=lw,label="TT-HEOM")	
	ax.plot(datN[2:end,1]*fac_t_E, (datN[2:end,2]-datN[2:end,3])*fac_c_E, color="red", linewidth=lw,label="Exact "*L"$t_1=0\,\mathrm{meV}$",linestyle="dashed")
	if i == 1
		ax.legend()
	end
	if i < length(Vs)
		ax.set_xticklabels([])
	else
		ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
	end
	ax.annotate(L"$V=$\, "*"$(round(V*0.2*0.04,digits=2))"*L"$\, \mathrm{V}$",(annx,anny),xycoords="axes fraction",fontsize=fs)
	ax.set_ylabel(L"I\,[\mu \mathrm{A}]")
	ax.set_ylim(yl[i][1],yl[i][2])
end
plt.savefig("Peierls_dyn_TTHEOM_Ehrenfest_noPh_VVar.pdf")


	

