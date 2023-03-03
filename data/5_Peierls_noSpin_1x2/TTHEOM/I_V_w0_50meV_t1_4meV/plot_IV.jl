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

V = [0,10,12,14]
Vlat = [L"$V=0\,\mathrm{V}$", 
	L"$V=0.4\,\mathrm{V}$",
	L"$V=0.48\,\mathrm{V}$",
	L"$V=0.56\,\mathrm{V}$",
	L"$V=0.8\,\mathrm{V}$"]
w0 = 0.05
t1 = 0.004 

ttIV = readdlm("I_V.dat")
mfIV = readdlm("../../MFHEOM/eps_-240_w0_50_t1_4/I_V.dat")
t0 = 0.04
ax.plot(ttIV[:,1]*t0,ttIV[:,2],color="k",label="TTHEOM",marker="x",linewidth=lw-2)
ax.plot(mfIV[:,1]*t0,mfIV[:,2]*2,color="red",label="MFHEOM",linewidth=lw-2,marker="o")

ax.legend(loc="center right")
ax.set_xlim(10*t0,18*t0)
ax.set_ylim(0.5,2.5)
ax.set_xlabel(L"$V\, [\mathrm{V}]$")
ax.set_ylabel(L"$I\,[\mu \mathrm{A}]$")
plt.savefig("Peierls_L2_w0_$(w0)_t1_$(t1)_I_V_comparison.pdf")
