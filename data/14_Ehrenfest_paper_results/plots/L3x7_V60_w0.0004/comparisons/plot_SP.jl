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

dat = "spin_selectivity_avg_n162.dat"

for i in 0:1
	data = readdlm("../runs/$(i)/$(dat)")
	datb = readdlm("../tsoPh0/$(i)/$(dat)")
	plt.plot(data[:,1]*fac_t,data[:,2],color="black",linewidth=lw)
	plt.plot(datb[:,1]*fac_t,datb[:,2],color="red",linewidth=lw)
end
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$S\,[\%]$")
plt.savefig("L3x7_V60_tsoPh_comp_SP.pdf")
