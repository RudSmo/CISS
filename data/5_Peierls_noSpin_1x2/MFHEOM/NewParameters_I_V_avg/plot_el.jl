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

H(X,eps,t0,t1,w) = [eps -t0-t1*sqrt(2*w)*X; -t0-t1*sqrt(2*w)*X eps]
eps = 0.0
t0  = 2.0
t1  = 0.2
w0  = 0.1


V = ["2.0","4.0","4.2","4.4","4.6","5.0"]
Vlat = [L"$2.0$",L"$4.0$",L"$4.2$",L"$4.4$",L"$4.6$",L"$5.0$"]
w = 0.1
colors = ["blue","green","red"]
for (v,vl,col) in zip(V,Vlat,colors)
	Vfold = Int(parse(Float64,v)*10)	
	dat = readdlm("V$(Vfold)/phonons_avg_n30.dat")
	t = dat[2:end,1]
	Xph = dat[2:end,2]

	el1 = []
	el2 = []

	for x in Xph
		Ham = H(x,eps,t0,t1,w0)
		ew = eigvals(Ham)
		push!(el1,ew[1])
		push!(el2,ew[2])
	end

	ax.plot(t*fac_t,el1,linewidth=lw,color=col,label=L"$V=\,$"*vl*" V")
	ax.plot(t*fac_t,el2,linewidth=lw,color=col,linestyle="dashed")
	muL = parse(Float64,v)/2
	ax.fill_between(t*fac_t, muL,-muL, color="lightskyblue",alpha=0.5)
#	ax.plot(t*fac_t,Xph,linewidth=lw,color=col)	
end

ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax.legend()
ax.set_ylabel(L"$E\,[\mathrm{eV}]$")
plt.savefig("Peierls_L2_NewParams_EL_VVar_avg.pdf")
