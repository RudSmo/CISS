using LinearAlgebra
using DelimitedFiles
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)
lw = 3
al = 0.4
fac_t = 1/(4.135668*1e3)
fac_c = 1e3 * 2.8494561
fac_p = 0.001

#plot = "electrons"
plot = "spin"
#phonons = false
phonons = true

N = 16
V = 0.48
tPh = 0.004
#tPh = 0.0
#tso = 0.001
tso = 0.0
#tsoPh = 0.0001
tsoPh = 0.0
#Tend = 12500
Tend = 37500
#Tend = 6250
ver = 2
#ver = 1
#P = 0.0
P = 0.5

dat_el = readdlm("L16_tPh0.004_SOC0.0_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.05_2_electrons.dat")
dat_ph = readdlm("L16_tPh0.004_SOC0.0_SOCPh0.0_V0.48_w00.002_P0.5_X0_-10_dt0.05_2_phonons.dat")

fig = plt.figure(figsize=(10,6))
px = 0.095
py = 0.55
sx = 0.8
sy = 0.4
dx = 0.02
dy = 0.04

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px,py-sy-dy,sx,sy])

ax1.plot(dat_el[2:end,1]*fac_t,0.5*(dat_el[2:end,10]-dat_el[2:end,11]),color="k",linewidth=lw)
ax1.set_ylabel(L"$s_5^z$")
ax1.set_ylim(-0.08,0.08)
ax1.set_xticklabels([])
ax1.annotate("(a)",(0.03,0.89),xycoords="axes fraction",fontsize=fs)
ax11 = ax1.twinx()


X5 = dat_ph[2:end,10]
X6 = dat_ph[2:end,12]
X7 = dat_ph[2:end,14]

ax11.plot(dat_ph[2:end,1]*fac_t,X6,label=L"$X_5$",linewidth=lw)
ax11.set_ylabel(L"$X_5$",color="blue")


ax2.plot(dat_el[2:end,1]*fac_t,0.5*(dat_el[2:end,12]-dat_el[2:end,13]),color="k",linewidth=lw)
ax2.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax2.set_ylabel(L"$s_6^z$")
ax2.set_ylim(-0.08,0.08)
ax2.annotate("(b)",(0.03,0.89),xycoords="axes fraction",fontsize=fs)
ax22 = ax2.twinx()
ax22.plot(dat_ph[2:end,1]*fac_t,X7,label=L"$X_6$",linewidth=lw)
ax22.set_ylabel(L"$X_6$",color="blue")


plt.savefig("L$(N)_site_6_el_ph.pdf")


