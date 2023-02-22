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

px = 0.09
py = 0.1
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

cur = readdlm("Peierls_L2_currents.dat")
ax.plot(cur[2:end,1],cur[2:end,2]*fac_c,color="k",linewidth=lw)
ax.set_xlabel(L"$V\,[V]$")
ax.set_ylabel(L"$I\,[\mu A]$")
plt.savefig("Peierls_L2_noPh_IV.pdf")

