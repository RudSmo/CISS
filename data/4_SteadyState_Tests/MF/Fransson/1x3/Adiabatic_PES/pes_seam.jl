using LinearAlgebra
using PyPlot
using3D()
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 1
num_ions = 3
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.1*t0
tsoPh = 0.1*tso
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5*HamPar.ws[i]^2 * HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.ws[i]^2 * HamPar.X[i]^2
        end
        return H
end


fig = plt.figure(figsize=(12,6))
px = 0.055
py = 0.1
sx = 0.415
sy = 0.88
dx = 0.015
dy = 0.04

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
#ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
axc = plt.axes([px+2*sx+dx+0.005,py,0.0125,sy])

Xrng = -250:1:150
nx = length(Array(Xrng))
colors = ["blue","red","green","gray","black","orange"]

fac_d = 1e3
DPES1 = zeros(nx,nx)
DPES2 = zeros(nx,nx)
#DPES3 = zeros(nx,nx)
for (i,X1) in enumerate(Xrng)
	for (j,X2) in enumerate(Xrng)
		HamPar = HamiltonianParams( num_laps,
					    num_ions,
					    N,
					    model,
					    modelStr,
					    E,
					    a,
					    c,
					    t0,
					    tPh,
					    tso,
					    tsoPh,
					    gPh,
					    ws,
					    [X1,X2,0.0],
					    pnIn,
					    tL,
					    tR)


		Hel = model(HamPar)
		Hph = PhononHamiltonian(HamPar)
		H = Hel .+ Hph
		ew = eigvals(H)
		d12 = abs(ew[3]-ew[1])
		d23 = abs(ew[5]-ew[3])
#		d13 = abs(ew[5]-ew[1])
		DPES1[i,j] = d12
		DPES2[i,j] = d23
#		DPES3[i,j] = d13
	end
end

cmaps = "inferno_r"
cp = plt.cm.inferno_r
vmn = 0.0*fac_d #minimum(DPES1)
vmx = 0.015*fac_d #maximum(DPES1)
norm = plt.matplotlib.colors.Normalize(vmin=vmn, vmax=vmx)
cmap = plt.cm.ScalarMappable(cmap=cp,norm=norm)

ax1.imshow(DPES1*fac_d,aspect="auto",extent=[Xrng[1],Xrng[end],Xrng[end],Xrng[1]],vmin=vmn,vmax=vmx,cmap=cmaps,interpolation="none")
ax1.set_ylabel(L"$X_1$",labelpad=-11)
ax1.set_xlabel(L"$X_2$")
ax1.annotate("(a)",(0.05,0.95),xycoords="axes fraction",fontsize=fs,color="white")
ax1.annotate(L"$i=1,j=3$",(0.05,0.025),xycoords="axes fraction",fontsize=fs,color="white")

ax2.imshow(DPES2*fac_d,aspect="auto",extent=[Xrng[1],Xrng[end],Xrng[end],Xrng[1]],vmin=vmn,vmax=vmx,cmap=cmaps,interpolation="none")
ax2.set_yticklabels([])
ax2.set_xlabel(L"$X_2$")
ax2.annotate("(b)",(0.05,0.95),xycoords="axes fraction",fontsize=fs,color="white")
ax2.annotate(L"$i=3,j=5$",(0.05,0.025),xycoords="axes fraction",fontsize=fs,color="white")

#ax3.imshow(DPES3*fac_d,aspect="auto",extent=[Xrng[1],Xrng[end],Xrng[end],Xrng[1]],vmin=vmn,vmax=vmx,cmap=cmaps,interpolation="none")
#ax3.set_yticklabels([])
#ax3.set_xlabel(L"$X_2$")
#ax3.annotate("(c)",(0.05,0.95),xycoords="axes fraction",fontsize=fs,color="white")
#ax3.annotate(L"$i=1,j=3$",(0.05,0.025),xycoords="axes fraction",fontsize=fs,color="white")

cbar = fig.colorbar(cmap, ax=ax2, cax = axc, aspect=90,pad=.001)
cbar.set_label(L"$\vert\varepsilon_i-\varepsilon_j\vert\,[\mathrm{meV}]$",labelpad=5)
cbar.ax.tick_params(labelsize=fs)

plt.savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh)_SOCPh$(tsoPh)_multiState_difference.pdf")
