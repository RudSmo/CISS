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
py = 0.55
sx = 0.43
sy = 0.44
dy = 0.03
dx = 0.02

fig = plt.figure(figsize=(12,8))
ax1 = plt.axes([px,py,sx,sy])
tax1 = ax1.twinx()
ax2 = plt.axes([px+sx+dx,py,sx,sy])
tax2 = ax2.twinx()
ax3 = plt.axes([px,py-sy-dy,sx,sy])
tax3 = ax3.twinx()
ax4 = plt.axes([px+sx+dx,py-sy-dy,sx,sy])
tax4 = ax4.twinx()


H(X,eps,t0,t1,w) = [eps -t0-t1*sqrt(2*w)*X; -t0-t1*sqrt(2*w)*X eps]
eps = 0.0
t0  = 2.0
t1  = 0.2
w0  = 0.1

params = [("4.0","0.0"),("20.0","0.0"),("2.0","-1.0"),("2.0","-2.0")]
folders = ["i","ii","iii","iv"]
axes = [ax1,ax2,ax3,ax4]
taxes = [tax1,tax2,tax3,tax4]
w = 0.1
colors = ["blue","green","red","orange"]
for (ax,tax,par,fold,col) in zip(axes,taxes,params,folders,colors)
	lab = "("*fold*")"
	dat = readdlm("$(fold)/phonons_avg_n30.dat")
	t = dat[2:end,1]
	Xph = dat[2:end,2]

	el1 = []
	el2 = []

	for x in Xph
		Ham = H(x,parse(Float64,par[2]),t0,t1,w0)
		ew = eigvals(Ham)
		push!(el1,ew[1])
		push!(el2,ew[2])
	end

	ax.plot(t*fac_t,el1,linewidth=lw,color=col)
	ax.plot(t*fac_t,el2,linewidth=lw,color=col,linestyle="dashed")
	ax.annotate(lab,(0.1,0.9),xycoords="axes fraction",fontsize=fs)
	ax.set_ylim(-5,5)
	tax.plot(t*fac_t,Xph,linewidth=lw-1,color="purple")
	tax.set_ylim(-5,5)	
	muL = parse(Float64,par[1])/2
	ax.fill_between(t*fac_t, muL,-muL, color="lightskyblue",alpha=0.5)
end

for ax in [tax2,tax4]
	ax.set_ylabel(L"X",color="purple")
end
	
for ax in [ax3,ax4]
	ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
end
for (ax,tax) in zip([ax1,ax3],[tax1,tax3])
	ax.set_ylabel(L"$E\,[\mathrm{eV}]$")
	tax.set_yticklabels([])
end
for ax in [ax1,ax2]
	ax.set_xticklabels([])
end
for ax in [ax2,ax4]
        ax.set_yticklabels([])
end
plt.savefig("Peierls_L2_NewParams_Cases_EL_ParVar_Ph_avg.pdf")
