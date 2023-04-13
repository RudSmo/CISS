using LinearAlgebra
using PyPlot
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 2

fac_t = 1/(4.135668)


rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.87
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

dat_t0 = readdlm("V_t0/L50_SOC0.0_V0.04_P0.0_noEps_dt0.05_electrons_es.dat")
dat_3t0_4 = readdlm("V_3t0_4/L50_SOC0.0_V0.03_P0.0_noEps_dt0.05_electrons_es.dat")
dat_4t0_3 = readdlm("V_4t0_3/L50_SOC0.0_V0.05333333333333334_P0.0_noEps_dt0.05_electrons_es.dat")

ax.set_xlabel(L"$\varepsilon\,[\mathrm{meV}]$")
ax.set_ylabel(L"$n(\varepsilon)$")
Erng = -40:80/(length(dat_t0[2,2:2:end])-1):40
ax.plot(Erng,(dat_3t0_4[2,2:2:end]+dat_3t0_4[2,3:2:end])/2,label="initial state",linewidth=lw)
ax.plot(Erng,(dat_3t0_4[end,2:2:end]+dat_3t0_4[end,3:2:end])/2,label=L"$V=30\,\mathrm{mV}$",linewidth=lw)
ax.plot(Erng,(dat_t0[end,2:2:end]+dat_t0[end,3:2:end])/2,label=L"$V=40\,\mathrm{mV}$",linewidth=lw)
ax.plot(Erng,(dat_4t0_3[end,2:2:end]+dat_4t0_3[end,3:2:end])/2,label=L"$V=53.33\,\mathrm{mV}$",linewidth=lw)
ax.legend(loc="upper right")
plt.savefig("Momentumdist_L50_TB_VVar_T100K.pdf")


