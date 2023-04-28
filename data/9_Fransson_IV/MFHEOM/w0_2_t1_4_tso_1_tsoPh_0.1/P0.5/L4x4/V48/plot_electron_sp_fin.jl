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

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,7))
ax = plt.axes([px,py,sx,sy])
num_sites = 4
num_laps = 4
N = num_laps * num_sites
V = 0.04 * 0.2 * 48
P = 0.5
V = 48 * 0.2 * 0.04
dat = readdlm("electrons_avg_n101.dat")[end,2:end]
datno = readdlm("noPh/L16_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V0.38400000000000006_w00.002_P0.5_noEps_Np15_dt0.05_1_electrons.dat")[end,2:end]
dat_tsoPh_0 = readdlm("../../../../w0_2_t1_4_tso_1_tsoPh_0/P0.5/V48/electrons_avg_n51.dat")[end,2:end]
dat_tsoPh_1 = readdlm("../../../../w0_2_t1_4_tso_1_tsoPh_1/P0.5/V48/electrons_avg_n51.dat")[end,2:end]
 
ax.plot(datno[1:2:end]-datno[2:2:end],label=L"$t_1= \lambda_1=0\, \mathrm{meV}$",linewidth=lw,color="red")
ax.plot(dat[1:2:end]-dat[2:2:end],label=L"$t_1=4\,\mathrm{meV},\, \lambda_1=0.1\, \mathrm{meV}$",linewidth=lw,color="k")
ax.plot(dat_tsoPh_0[1:2:end]-dat_tsoPh_0[2:2:end],label=L"$t_1=4\,\mathrm{meV},\, \lambda_1=0\, \mathrm{meV}$",linewidth=lw,color="blue")
ax.plot(dat_tsoPh_1[1:2:end]-dat_tsoPh_1[2:2:end],label=L"$t_1=4\,\mathrm{meV},\, \lambda_1=1\, \mathrm{meV}$",linewidth=lw,color="green")

ax.set_xlabel("lattice site")
ax.set_ylabel(L"$s_i^z$")
ax.legend()


plt.savefig("L$(num_laps)x$(num_sites)_P$(P)_V$(V)_electron_sp_comp.png",dpi=200)
