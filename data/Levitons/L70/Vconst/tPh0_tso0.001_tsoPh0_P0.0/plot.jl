using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3

fac_t = 1/(4.135668)


rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.81
sxc = 0.01
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+sx+dx,py,sxc,sy])


for obs in ["local_nd","local_sz","momentum_occ_diff"]
	if obs == "local_nd" || obs == "local_sz"
		datname = readdir()[2] 
	elseif obs == "momentum_occ_diff"
		datname = readdir()[3]
	end
	splitname = split(datname,"_")
	tPh = splitname[3]
	SOC = splitname[6][4:end]
	SOCPh = splitname[7][6:end]
	P = splitname[11][2:end]

	dat = readdlm(datname)

	if obs == "local_nd"
		OBS = (dat[2:end,2:2:end]+dat[2:end,3:2:end])/2
		ext_r = 100
	elseif obs == "local_sz"
		OBS = 1e3*(dat[2:end,2:2:end]-dat[2:end,3:2:end])
		ext_r = 100
	elseif obs == "momentum_occ_diff"
		OBS = zeros(size(dat[2:end,2:end]))
		for i in 2:length(OBS[2:end,1])
			OBS[i-1,:] = dat[i,2:end] - dat[2,2:end]
		end
		ext_r = length(dat[2,2:end])
	end

	cp = plt.cm.jet
	vmn = minimum(OBS)
	vmx = maximum(OBS)
	norm = plt.matplotlib.colors.Normalize(vmin=vmn, vmax=vmx)
	cmap = plt.cm.ScalarMappable(cmap=cp,norm=norm)

	ax.imshow(OBS,aspect="auto",extent=[1,ext_r,fac_t*dat[end,1],fac_t*dat[2,1]],cmap="jet")
	im = ax.set_ylabel(L"$t\,[\mathrm{fs}]$")
	cbar = fig.colorbar(cmap,ax=ax,cax=cax,pad=0.001)
	if obs == "local_nd"
		cbar.set_label(L"$n$")
		ax.set_xlabel("site")
	elseif obs == "local_sz"
		cbar.set_label(L"$s^z\,[\times 10^{-3}]$",labelpad=-3)
		ax.set_xlabel("site")
	elseif obs == "momentum_occ_diff"
		cbar.set_label(L"$n(\varepsilon_k,t)-n(\varepsilon_k,t_0)$")
		ax.set_xlabel("state "*L"$\varepsilon_k$")
	end
	
	plt.savefig("Leviton_tPh_$(tPh)_SOC_$(SOC)_SOCPh_$(SOCPh)_P$(P)_$(obs).pdf")
end
