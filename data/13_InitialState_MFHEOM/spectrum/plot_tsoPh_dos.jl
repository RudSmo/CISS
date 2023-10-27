using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
lw = 3
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.11
py = 0.12
sx = 0.87
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(7,5))
ax = plt.axes([px,py,sx,sy])

#  2    3      4       5     6    7
# [0.0,0.0001,0.0005,0.001,0.005,0.01]
dat = readdlm("L2x4_tPh_[0.0, 0.0]_tsoPh_var_avg_dos.dat")
ax.plot(dat[:,1],dat[:,3] .+ 1000,label=L"$t_1 = 0, \lambda_1=0.1\, \mathrm{meV}$",linewidth=lw,color="k")
ax.plot(dat[:,1],dat[:,5] .+ 700,label=L"$\lambda_1=1\, \mathrm{meV}$",linewidth=lw,color="blue")
ax.plot(dat[:,1],dat[:,6] .+ 200,label=L"$\lambda_1=5\, \mathrm{meV}$",linewidth=lw,color="red")
ax.plot(dat[:,1],dat[:,7],label=L"$\lambda_1=10\, \mathrm{meV}$",linewidth=lw,color="green")
ax.legend(fontsize=fs-2,loc="lower left",framealpha=0.5)
ax.set_ylim(0,1500)
ax.set_xlabel(L"$E\,[\mathrm{meV}]$")
ax.set_ylabel("DOS")
plt.savefig("L2x4_tPh_0.0_tsoPhVar_dos.pdf")

