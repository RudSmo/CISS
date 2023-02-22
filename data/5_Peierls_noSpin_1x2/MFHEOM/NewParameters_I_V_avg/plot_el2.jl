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

px = 0.06
py = 0.095
sx = 0.28
sy = 0.88
dy = 0.03
dx = 0.02

fig = plt.figure(figsize=(12,6))
ax1 = plt.axes([px,py,sx,sy])
tax1 = ax1.twinx()
ax2 = plt.axes([px+sx+dx,py,sx,sy])
tax2 = ax2.twinx()
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
tax3 = ax3.twinx()


H(X,eps,t0,t1,w) = [eps -t0-t1*sqrt(2*w)*X; -t0-t1*sqrt(2*w)*X eps]
eps = 0.0
t0  = 2.0
t1  = 0.2
w0  = 0.1

V = ["0.0","4.0","7.0"]
Vlat = [L"$0.0$",L"$4.0$",L"$7.0$"]
axes = [ax1,ax2,ax3]
taxes = [tax1,tax2,tax3]
w = 0.1
colors = ["blue","green","red"]
for (ax,tax,v,vl,col) in zip(axes,taxes,V,Vlat,colors)
	dat = readdlm("L2_gPh0.0_tPh-0.2_V$(v)_w00.1_E0.0_dt0.0005_1_phonons.dat")
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

	ax.plot(t*fac_t,el1,linewidth=lw,color=col)
	ax.annotate(L"$V=\,$"*vl*" V",(0.1,0.9),xycoords="axes fraction",fontsize=fs)
	ax.plot(t*fac_t,el2,linewidth=lw,color=col)
	ax.set_ylim(-5,5)
	ax.set_xlim(0,260)
	tax.plot(t*fac_t,Xph,linewidth=lw-1,color="purple")
	tax.set_ylim(-10,10)	
	muL = parse(Float64,v)/2
	ax.fill_between(t*fac_t, muL,-muL, color="lightskyblue",alpha=0.5)
end

tax3.set_ylabel(L"$X \, [\mathrm{eV}^{-1/2}]$",color="purple",labelpad=-3)
	
for ax in [ax1,ax2,ax3]
	ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
end
ax1.set_ylabel(L"$E\,[\mathrm{eV}]$")

for (tax,ax) in zip([tax1,tax2],[ax2,ax3])
	tax.set_yticklabels([])
	ax.set_yticklabels([])
end
plt.savefig("Peierls_L2_NewParams_EL_ParVar_Ph.pdf")
