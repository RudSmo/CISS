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
py = 0.54
sx = 0.27
sy = 0.44
dy = 0.03
dx = 0.02

fig = plt.figure(figsize=(12,8))
ax1 = plt.axes([px,py,sx,sy])
tax1 = ax1.twinx()
ax2 = plt.axes([px+sx+dx,py,sx,sy])
tax2 = ax2.twinx()
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
tax3 = ax3.twinx()
ax4 = plt.axes([px,py-sy-dy,sx,sy])
tax4 = ax4.twinx()
ax5 = plt.axes([px+sx+dx,py-sy-dy,sx,sy])
tax5 = ax5.twinx()
ax6 = plt.axes([px+2*(sx+dx),py-sy-dy,sx,sy])
tax6 = ax6.twinx()



H(X,eps,t0,t1,w) = [eps -t0-t1*sqrt(2*w)*X; -t0-t1*sqrt(2*w)*X eps]
eps = 0.0
t0  = 2.0
t1  = 0.2
w0  = 0.1

params = [2.0,4.0,4.2,4.4,4.6,5.0]
labels = [L"$V=2.0\,\mathrm{V}$",L"$V=4.0\,\mathrm{V}$",L"$V=4.2\,\mathrm{V}$",L"$V=4.4\,\mathrm{V}$",L"$V=4.6\,\mathrm{V}$",L"$V=5.0\,\mathrm{V}$"]
colors = ["blue","green","red","orange","gray","magenta"]
axes = [ax1,ax2,ax3,ax4,ax5,ax6]
taxes = [tax1,tax2,tax3,tax4,tax5,tax6]
w = 0.1
for (ax,tax,par,lab,col) in zip(axes,taxes,params,labels,colors)
	Vfold = Int(par*10)
	datPh = readdlm("V$(Vfold)/phonons_avg_n30.dat")
	datES = readdlm("V$(Vfold)/electrons_es_avg_n30.dat")

	t = datPh[2:end,1]
	Xph = datPh[2:end,2]
	
	Dif = datES[2:end,2] - datES[2:end,3]

	ax.plot(t*fac_t,t1*sqrt(2/w0^3) * Dif,linewidth=lw,color=col)
	ax.annotate(lab,(0.025,0.92),xycoords="axes fraction",fontsize=fs)
	if ax == ax1 || ax == ax2 || ax == ax3
		ax.set_ylim(-0.1,10)
		tax.set_ylim(-0.1,10)
	else
		ax.set_ylim(-0.2,1.4)
		tax.set_ylim(-0.2,1.4)
	end
	tax.plot(t*fac_t,Xph,linewidth=lw-1,color="purple")
end

for ax in [tax3,tax6]
	ax.set_ylabel(L"$X$",color="purple")
end
	
for ax in [ax4,ax5,ax6]
	ax.set_xlabel(L"$t\,[\mathrm{fs}]$")
end
for ax in [ax1,ax4]
	
	if ax == ax4
		lp = -2
	else
		lp = 7
	end
	ax.set_ylabel(L"$t_1\sqrt{\frac{2}{\omega^3}}(\vert c_{-}\vert^{2} - \vert c_{+}\vert^{2})$",labelpad=lp)
end
for tax in [tax1,tax2,tax4,tax5]
	tax.set_yticklabels([])
end
for ax in [ax1,ax2,ax3]
	ax.set_xticklabels([])
end
for ax in [ax2,ax3,ax5,ax6]
        ax.set_yticklabels([])
end
plt.savefig("Peierls_L2_NewParams_EL_ES_ParVar_Ph_avg.pdf")
