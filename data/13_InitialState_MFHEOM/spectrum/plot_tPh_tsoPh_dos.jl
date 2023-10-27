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

dat = readdlm("L2x4_tPh_var_tsoPh_var_avg_dos.dat")
ax.plot(dat[:,1],dat[:,2] .+ 900,label=L"$\lambda_1/\lambda_0 = t_1/t_0 = 0$",linewidth=lw,color="k")
ax.plot(dat[:,1],dat[:,3] .+ 700,label=L"$0.1$",linewidth=lw,color="blue")
ax.plot(dat[:,1],dat[:,4] .+ 400,label=L"$0.2$",linewidth=lw,color="red")
ax.plot(dat[:,1],dat[:,5],label=L"$1$",linewidth=lw,color="green")
ax.legend(fontsize=fs-2,loc="lower left",framealpha=0.5)
ax.set_ylim(0,1100)
ax.set_xlabel(L"$E\,[\mathrm{meV}]$")
ax.set_ylabel("DOS")
plt.savefig("L2x4_tsoPhVar_tPhVar_dos.pdf")

