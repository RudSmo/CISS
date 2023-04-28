using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.14
sx = 0.87
sy = 0.83
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,5))
ax = plt.axes([px,py,sx,sy])
num_sites = 4
num_laps = 4
N = num_laps * num_sites
V = 0.04 * 0.2 * 48
P = 0.5
Vs = ["V48","V55","V65"]
labels = [L"$V=\,0.384\mathrm{V}$",L"$V=\,0.44\mathrm{V}$",L"$V=\,0.52\mathrm{V}$"]
colors = ["black","blue","green"]
colors2 = ["brown","teal","magenta"]
for (v,lab,col,col2) in zip(Vs,labels,colors,colors2)
        dat = readdlm("$(v)/currents_avg_n101.dat")
	Id = (dat[:,2]-dat[:,3])*fac_c
	dat_std = readdlm("$(v)/currents_std_n101.dat")
	cur_std = sqrt.((dat_std[:,2]*fac_c) .^ 2 .+ (dat_std[:,3]*fac_c) .^ 2)
	cur_serr = cur_std/sqrt(101)

	curfile = fnfilter(readdir("$(v)/noPh/"),"*currents*")[1]
	datNo = readdlm("$(v)/noPh/$(curfile)")

	
        ax.plot(dat[:,1]*fac_t*1e-3,Id,linewidth=lw,label=lab,color=col)	
	ax.fill_between(dat[:,1]*fac_t*1e-3,Id,Id+cur_serr,color=col2,alpha=0.5)
	ax.fill_between(dat[:,1]*fac_t*1e-3,Id,Id-cur_serr,color=col2,alpha=0.5)
#	ax.plot(datNo[2:end,1]*fac_t*1e-3,(datNo[2:end,2]-datNo[2:end,3])*fac_c,linewidth=lw-1,linestyle="dashed",color=col)
#	ax.plot(datNo[end-100,1]*fac_t*1e-3:1:31,ones(length(datNo[end-100,1]*fac_t*1e-3:1:31))*(datNo[end,2]-datNo[end,3])*fac_c,linewidth=lw-1,linestyle="dashed",color=col)
end
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
ax.set_xlim(-0.1,30.1)
ax.set_ylim(0,1.4)
ax.legend()


plt.savefig("L$(num_laps)x$(num_sites)_P$(P)_VVar_current_101_noComp.png",dpi=200)
