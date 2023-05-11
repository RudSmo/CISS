using PyPlot
using LinearAlgebra
include("../../../src/MeanField_v3/classical_phonons.jl")
using .ClassicalPhonons

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

px = 0.09
py = 0.58
sx1 = 0.37
sx = 0.4
sy = 0.4
dy = 0.1
dx = 0.08
dx2 = 0.02

fig = plt.figure(figsize=(9,7))
ax1 = plt.axes([px,py,sx1,sy])
ax2 = plt.axes([px+sx1+dx-0.03,py,sx1,sy])
cax2 = plt.axes([px+2*sx1+dx-0.02,py,0.013,sy])
ax3 = plt.axes([px,py-sy-dy,sx-dx2,sy])
ax4 = plt.axes([px+sx+dx+dx2,py-sy-dy,sx-dx2,sy])


function Wigner(X,P,beta,w)
	nX = length(X)
	nP = length(P)
	G = zeros(nP,nX)
	a = tanh(beta*w/2)
	for i in 1:nX
		for j in 1:nP
			G[j,i] = 2*a*exp(-a*(X[i]^2/2+P[j]^2/2))
		end
	end
	return G
end

beta = 40
w = 0.002
ed = 20
X = -ed:0.01:ed
P = -ed:0.01:ed

# analytical distribution
WA_full = Wigner(X,P,beta,w)
WA_X = reshape(Wigner(X,[0.0],beta,w),length(X))
WA_P = reshape(Wigner([0.0],P,beta,w),length(P))

# sampled Wigner distribution
WS_X,WS_P = SampleWigner_Multi(beta,w,100000)



cp = plt.cm.inferno
cmap = plt.cm.ScalarMappable(cmap=cp)
vmn = 0.0 #minimum(WA_full)
vmx = maximum(WA_full)

im = ax1.imshow(WA_full,aspect="auto",vmin=vmn,vmax=vmx,extent=[-ed,ed,-ed,ed],cmap="inferno")
ax1.set_ylabel(L"$P$")
ax1.set_xlabel(L"$X$")
ax1.annotate("analytical",(0.05,0.1),xycoords="axes fraction",color="white",fontsize=fs)
annx = 0.03
anny = 0.91
ax1.annotate("(a)",(annx,anny),xycoords="axes fraction",color="white",fontsize=fs)


ax2.hist2d(WS_X,WS_P,density=true,bins=40,cmap="inferno")

ax2.set_yticklabels([])
ax2.annotate("sampled",(0.05,0.1),xycoords="axes fraction",color="white",fontsize=fs)
ax2.set_xlim(-ed,ed)
ax2.set_xlabel(L"$X$")
ax2.annotate("(b)",(annx,anny),xycoords="axes fraction",color="white",fontsize=fs)
cbar = fig.colorbar(im,ax=ax2,cax = cax2,aspect=90,pad=.001)
cbar.set_label(L"$\rho_W^\mathrm{G}(X,P)$")
cbar.ax.tick_params(labelsize=fs)
cbar.set_ticks(Array(0.0:0.01:0.08))
cbar.set_ticklabels([string(i) for i in 0.0:0.01:0.08])

ax3.hist(WS_X,density=true,bins=100)
ax3.plot(X,WA_X,linewidth=lw,color="orange",alpha=0.7)
ax3.set_ylim(-0.003,0.083)
ax3.set_xlim(-ed,ed)
ax3.set_ylabel(L"$\rho_W^\mathrm{G}(X,0)$")
ax3.set_xlabel(L"$X$")
ax3.annotate("(c)",(annx,anny),xycoords="axes fraction",fontsize=fs)

ax4.hist(WS_P,density=true,bins=100)
ax4.plot(P,WA_P,linewidth=lw,color="orange",alpha=0.7)
ax4.set_ylim(-0.003,0.083)
ax4.set_xlim(-ed,ed)
ax4.set_ylabel(L"$\rho_W^\mathrm{G}(0,P)$")
ax4.set_xlabel(L"$P$")
ax4.annotate("(d)",(annx,anny),xycoords="axes fraction",fontsize=fs)
plt.savefig("Comparison_sampled_Wigner.pdf")

