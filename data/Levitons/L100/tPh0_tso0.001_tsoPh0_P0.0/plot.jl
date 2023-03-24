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

datname = readdir()[1] 
splitname = split(datname,"_")

tPh = splitname[3]
SOC = splitname[6][4:end]
SOCPh = splitname[7][6:end]
P = splitname[11][2:end]

for obs in ["local_nd","local_sz"]
	dat = readdlm(datname)

	if obs == "local_nd"
		OBS = (dat[2:end,2:2:end]+dat[2:end,3:2:end])/2
	elseif obs == "local_sz"
		OBS = 1e3*(dat[2:end,2:2:end]-dat[2:end,3:2:end])
	end

	cp = plt.cm.jet
	vmn = minimum(OBS)
	vmx = maximum(OBS)
	norm = plt.matplotlib.colors.Normalize(vmin=vmn, vmax=vmx)
	cmap = plt.cm.ScalarMappable(cmap=cp,norm=norm)

	ax.imshow(OBS,aspect="auto",extent=[1,100,fac_t*dat[end,1],fac_t*dat[2,1]],cmap="jet")
	ax.set_xlabel("site")
	im = ax.set_ylabel(L"$t\,[\mathrm{fs}]$")
	cbar = fig.colorbar(cmap,ax=ax,cax=cax,pad=0.001)
	if obs == "local_nd"
		cbar.set_label(L"$n$")
	elseif obs == "local_sz"
		cbar.set_label(L"$s^z\,[\times 10^{-3}]$",labelpad=-3)
	end
	plt.savefig("Leviton_tPh_$(tPh)_SOC_$(SOC)_SOCPh_$(SOCPh)_P$(P)_$(obs).pdf")
end
