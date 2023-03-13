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
fac_t = 1e-3/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

w = 0.002
t1 = 0.004
L = "L2x3" 
V = [48,55,65]
for v in V
	datP = readdlm("P0.5/$(L)/V$(v)/currents_avg_n20.dat")
	datM = readdlm("P-0.5/$(L)/V$(v)/currents_avg_n20.dat")
	CurP = datP[2:end,2] .- datP[2:end,3]
	CurM = datM[2:end,2] .- datM[2:end,3]
	Dcur = (CurP .- CurM) * 100 ./ (CurP .+ CurM)
	ax.plot(datP[2:end,1]*fac_t, Dcur, label = L"$V=\,$"*string(0.2*0.04*v))	
end
ax.legend()
ax.set_ylim(-0.01,0.01)
ax.set_xlabel(L"$t\, [\mathrm{ps}]$")
ax.set_ylabel(L"$M_R\, [\%]$")

plt.savefig("$(L)_w0_$(w)_t1_$(t1)_Magnetocurrent_Vvar.pdf")


