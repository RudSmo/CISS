using PyPlot
using DelimitedFiles


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])

datPn = "convergence_info_avg_n101_P.dat"
datMn = "convergence_info_avg_n101_M.dat"

for i in 0:4
	datP = readdlm("../$(i)/$(datPn)")
	datM = readdlm("../$(i)/$(datMn)")
	plt.plot(datP[:,1]*fac_t,datP[:,3]*0.002*0.2,linewidth=lw,color="red")
	plt.plot(datM[:,1]*fac_t,datM[:,3]*0.002*0.2,linewidth=lw,color="blue")
end
ax.set_ylabel(L"$E_\mathrm{sys}\,[\mathrm{eV}]$")
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
plt.savefig("L3x7_V60_Esys_PM.pdf")
