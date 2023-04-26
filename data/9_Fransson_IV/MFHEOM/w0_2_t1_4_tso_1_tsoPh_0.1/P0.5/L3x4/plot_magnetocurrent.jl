using LinearAlgebra
using PyPlot
using DelimitedFiles

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

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
num_sites = 4
num_laps = 3
N = num_laps * num_sites
Vs = ["V48","V55","V65"]
labels = [L"$V=\,0.384\mathrm{V}$",L"$V=\,0.44\mathrm{V}$",L"$V=\,0.52\mathrm{V}$"]
colors = ["blue","red","black"]
for (v,lab,col) in zip(Vs,labels,colors)
        datP = readdlm("$(v)/currents_avg_n51.dat")
        datM = readdlm("../../P-0.5/L$(num_laps)x$(num_sites)/$(v)/currents_avg_n51.dat")

	datstdP = readdlm("$(v)/currents_std_n51.dat") 
	stdPL   = datstdP[:,2]
	stdPR   = datstdP[:,3] 
	datstdM = readdlm("../../P-0.5/L$(num_laps)x$(num_sites)/$(v)/currents_std_n51.dat")
	stdML   = datstdM[:,2]
	stdMR   = datstdM[:,3]

	curT = datP[:,2]-datP[:,3] + datM[:,2]-datM[:,3]
	curD = datP[:,2]-datP[:,3] - (datM[:,2]-datM[:,3])

	dPL = (curT .- curD) ./ (curT.^2)
        dPR = -dPL
        dML = -(curT .+ curD) ./ (curT.^2)
        dMR = -dML
        mr_serr = 100/sqrt(101) * sqrt.( (dPL .* stdPL).^2 .+ (dPR .* stdPR).^2 .+ (dML .* stdML).^2 + (dMR .* stdMR) .^ 2 )
	
	MR = 100 * curD ./ curT
        ax.plot(datP[:,1]*fac_t*1e-3,MR,color=col,linewidth=lw,label=lab)
	ax.fill_between(datP[:,1]*fac_t*1e-3,MR,MR+mr_serr,color=col,alpha=0.2)
        ax.fill_between(datP[:,1]*fac_t*1e-3,MR,MR-mr_serr,color=col,alpha=0.2)
	
end

ax.set_ylim(-25,25)
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$M_R\,[\%]$")
ax.legend()
plt.savefig("L$(num_laps)x$(num_sites)_VVar_magnetocurrent.pdf")
