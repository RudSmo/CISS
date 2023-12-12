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
py = 0.12
sx = 0.84
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,5))

t_E_0 = 1e3 * 0.2433946
t_E_1 = 1/(41.35668)
t_E_2 = 1/(0.002*41.35668)
I_E_1 = 1e3 * 0.2433946
I_E_2 = 0.002 * 1e3 * 0.2433946
t_T = 1e-3/2.418884326505E-2
ylim = [[1.1,1.3],[0.8,1.3],[0.3,0.75]]
ax = plt.axes([px,py,sx,sy])
files = ["tamen_V0.48_conv","tamen_V0.48_mvib30","tamen_V0.48_mvib50","tamen_V0.48_mvib70"]
labels = [L"$n_{vib}=20,\ \mathrm{res\_damp} = 100$", L"$n_{vib}=30,\ \mathrm{res\_damp} = 10$", L"$n_{vib}=50,\ \mathrm{res\_damp} = 10$", L"$n_{vib}=70,\ \mathrm{res\_damp} = 10$"]
for (f,lab) in zip(files,labels)
	dat = readdlm("$(f)/nvib.txt")
	sorting = sortperm(dat[:,1])
	dat = dat[sorting,:]
	TAMEN = unique(dat,dims=1)
	ax.plot(TAMEN[:,1],TAMEN[:,2],label=lab,linewidth=lw)
end
ax.legend(loc="lower right")
ax.set_xlim(0,600)
#ax.set_ylim(1.1,1.25)
ax.set_ylabel(L"$\langle a^\dagger a\rangle$")
ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
plt.savefig("Peierls_L2_spin_V60_mvibVar_nvib.pdf")



	


