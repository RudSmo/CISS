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

px = 0.07
py = 0.135
sx = 0.91
sy = 0.835
dx = 0.05

fig = plt.figure(figsize=(10,4))
ax = plt.axes([px,py,sx,sy])

V = [0,10,12,14,20]
Vlat = [L"$V=0.0\, \mathrm{V}$",
	L"$V=0.4\, \mathrm{V}$",
	L"$V=0.48\, \mathrm{V}$",
	L"$V=0.56\, \mathrm{V}$",
	L"$V=0.8\,\mathrm{V}$"]
w = 0.1 
colors = ["red",
	  "blue",
	  "green",
	  "orange",
	  "gray"]
for (v,vl,col) in zip(V,Vlat,colors)
	datQ = readdlm("V$(v)_t1_0.004_w0_0.002_nv50.dat")
	ax.plot(datQ[end,4:end],linewidth=lw,color=col,label=vl)
end

ax.legend(loc="upper right")
ax.set_xlim(0,20)
ax.set_ylabel(L"$\vert \langle n \vert\psi(t_\mathrm{f})\rangle\vert^2$")
ax.set_xlabel("Vibrational quantum number "*L"$n$")

plt.savefig("Peierls_L2_FranssonParams_PhQ_nv.pdf")
