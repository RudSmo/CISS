using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)
ws = 0.002
fac_X = 1.9732*1e-7 * 1e10 * ws
fac_P = 5.3442883*1e-28 * 1e-12 / (sqrt(ws) * 6.022*1e26 * 1e-10)
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

N = 5
dat = readdlm("V48/phonons_comb.dat")

ax = plt.figure().add_subplot(projection="3d")

t = 0:10:125000
for i in 1:N
	X = dat[end:-1:1,2*i]
	P = dat[end:-1:1,2*i+1]
	ax.plot(X,P,t)
end
plt.show()
	

