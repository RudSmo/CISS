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

px = 0.11
py = 0.55
sx = 0.86
sy = 0.435
dy = 0.04

fig = plt.figure(figsize=(10,8))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

datP = readdlm("P0.5/GF_DOS_N_8_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_ws_0.002_beta_40.0_tl_0.004_nS_10000_P_0.5_coupled_avg_100000.dat") 
datM = readdlm("P-0.5/GF_DOS_N_8_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_ws_0.002_beta_40.0_tl_0.004_nS_10000_P_-0.5_coupled_avg_100000.dat") 


Energy = datP[:,1]
dE = Energy[2]-Energy[1]
Dif = datP[:,2] .- datM[:,2] 
ax1.plot(Energy,Dif, color="k",linewidth=lw)
ax1.set_ylabel(L"$\mathrm{DOS}_\uparrow -\mathrm{DOS}_\downarrow $")
ax1.set_xticklabels([])
ax1.annotate("(a)",(0.02,0.91),xycoords="axes fraction",fontsize=fs)
ax1.set_xlim(-0.35,-0.15)

function Integrate(f,dx)
	f_out = []
	fint = 0.0
	N = length(f)
	for i in 1:N-1
		fint += 0.5 * dx * (f[i+1]+f[i])
		push!(f_out,fint)
	end
	return f_out
end


int = Integrate(Dif,dE)
ax2.plot(Energy[1:end-1],int,color="k",linewidth=lw)
ax2.set_ylabel(L"$\int_{-\infty}^\varepsilon \mathrm{d}\varepsilon' (\mathrm{DOS}_\uparrow -\mathrm{DOS}_\downarrow)$")
ax2.set_xlabel(L"$\varepsilon\,[\mathrm{eV}]$")
ax2.set_xlim(-0.35,-0.15)
ax2.annotate("(b)",(0.02,0.91),xycoords="axes fraction",fontsize=fs)
plt.savefig("L2x4_Fransson_diff_DOS_avg_1e5.pdf")
