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

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
axi = plt.axes([px+0.06,py+0.48,sx-0.56,sy-0.55])

datP = readdlm("P0.5/GF_DOS_N_8_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_ws_0.002_beta_40.0_tl_0.004_nS_1000_P_0.5_coupled_avg_11000.dat") 
datM = readdlm("P-0.5/GF_DOS_N_8_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_ws_0.002_beta_40.0_tl_0.004_nS_1000_P_-0.5_coupled_avg_11000.dat") 
datN = readdlm("P0.0/GF_DOS_N_8_eps_-0.24_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_ws_0.002_beta_40.0_tl_0.004_nS_1000_P_0.0_coupled_avg_11000.dat")
dat0 = readdlm("../tPh0.0_tso0.001_tsoPh0.0/GF_DOS_N_8_tso_0.001_tPhU_0.0_tPhD_0.0_tsoPh_0.0_beta_40.0_tl_0.004_nS_1_P_0.0_coupled.dat")

for axe in [ax,axi]
	axe.plot(datP[:,1],datP[:,2],color="red",linewidth=lw,label=L"$P=0.5$")
	axe.plot(datM[:,1],datM[:,2],color="blue",linewidth=lw,label=L"$P=-0.5$")
	axe.plot(datN[:,1],datN[:,2],color="k",linewidth=lw,label=L"$P=0.0$")
	axe.plot(dat0[:,1],dat0[:,2],color="orange",linewidth=lw-2,alpha=0.7,label=L"$P=0.0, \lambda_1 = t_1 = 0\,\mathrm{eV}$")
end

ax.set_xlabel(L"$\varepsilon\,[\mathrm{eV}]$")
ax.set_ylabel(L"$\mathrm{DOS}(\varepsilon)$")
ax.set_ylim(0,1000)
axi.set_xlim(-0.323,-0.312)
axi.set_ylim(120,160)
ax.legend()
plt.savefig("L2x4_Fransson_DOS_avg_1e4.pdf")
