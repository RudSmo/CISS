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

w = 0.01
t1 = 0.02


for i in 1:15
#	dat = readdlm("data$(i)/L2_tPhU_-0.04_tPhD_-0.04_SOC0.0_SOCPh0.0_V0.48_dt0.05_1_currents.dat")
#	ax.plot(dat[2:end,1]*fac_t,(dat[2:end,2]-dat[2:end,3])*fac_c,linewidth=lw-2,alpha=0.5)

#	dat = readdlm("data$(i)/L2_tPhU_-0.04_tPhD_-0.04_SOC0.0_SOCPh0.0_V0.48_dt0.05_1_phonons.dat")
#	ax.plot(dat[2:end,2],dat[2:end,3],linewidth=lw-2,alpha=0.5)
	
	dat = readdlm("data$(i)/L2_tPhU_-0.04_tPhD_-0.04_SOC0.0_SOCPh0.0_V0.48_dt0.05_1_electrons_es.dat")
	ax.plot(dat[2:end,1]*fac_t,(dat[2:end,2]+dat[2:end,3])/2,linewidth=lw-2,alpha=0.5)
	ax.plot(dat[2:end,1]*fac_t,(dat[2:end,4]+dat[2:end,5])/2,linewidth=lw-2,linestyle="dashed",alpha=0.5)
end

#ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
#ax.set_ylabel(L"$I \,[\mu\mathrm{A}]$")
#ax.set_ylim(0,2.)
#ax.set_xscale("log")
#plt.savefig("Fransson_L2_P0_non-unique_steady_state_currents.pdf")

#ax.set_xlabel(L"$X$")
#ax.set_ylabel(L"$P$")
#plt.savefig("Fransson_L2_P0_non-unique_steady_state_phonons.pdf")

ax.set_xlabel(L"$t\, [\mathrm{fs}]$")
ax.set_ylabel("Eigenstate population")
ax.set_xscale("log")
plt.savefig("Fransson_L2_P0_non-unique_steady_state_el_es.pdf")


