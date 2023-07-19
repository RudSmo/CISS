using LinearAlgebra
using PyPlot
using3D()
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
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


num_laps = 2
num_ions = 2
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.0 * T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.0 * 0.1*tso
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


fig = plt.figure(figsize=(10,6))
px = 0.09
py = 0.56
sx = 0.29
sy = 0.41
dx = 0.01
dy = 0.04

ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
ax4 = plt.axes([px,py-sy-dy,sx,sy])
ax5 = plt.axes([px+sx+dx,py-sy-dy,sx,sy])
ax6 = plt.axes([px+2*(sx+dx),py-sy-dy,sx,sy])

Xrng = -20:0.1:20
nx = length(Array(Xrng))
colors = ["blue","red","green","gray","black","orange"]
#X2s = [-400,-300,-200,-100,0,300]
X2s  = [-20,-10,0,10,15,20]
axes = [ax1,ax2,ax3,ax4,ax5,ax6]

for (X2,ax) in zip(X2s,axes)
	PES = []
	for (i,X1) in enumerate(Xrng)
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
					    [X1,X2,X2,0.0],
					    pnIn,
					    tL,
					    tR)


			Hel = model(HamPar)
			Hph = PhononHamiltonian(HamPar)
			H = Hel .+ Hph
			ew = eigvals(H)
			push!(PES,ew)
	end
	PES = [[PES[i][j] for i in 1:length(PES)] for j in 1:length(PES[1])]
	if X2 < 0
		ax.annotate(L"$X_2=-$"*string(abs(X2)),(0.05,0.05),xycoords="axes fraction",fontsize=fs)
	else
		ax.annotate(L"$X_2=\,$"*string(X2),(0.05,0.05),xycoords="axes fraction",fontsize=fs)
	end
	for (pes,col) in zip(PES,colors)
		ax.plot(Xrng,pes,linewidth=lw,color=col)
	end
	if ax == ax1 || ax == ax4
		ax.set_ylabel(L"$\varepsilon\, [\mathrm{eV}]$")
	else
		ax.set_yticklabels([])
	end
	if ax == ax1 || ax == ax2 || ax == ax3
		ax.set_xticklabels([])
	else
		ax.set_xlabel(L"$X_1$")
		#ax.set_xticklabels([L"$-500$",L"$-250$",L"$0$",L"$250$",L"$500$",L"$750$"])
	end
#	if ax == ax4
#		ax.set_xticklabels([-500,-250,0,250,500,750])
#	elseif ax == ax5
#		ax.set_xticklabels([-500,-250,0,250,500,750])
#	elseif ax == ax6
#		ax.set_xticklabels([-500,-250,0,250,500,750,1000])
#	end
	
#	ax.set_ylim(-0.4,0.2)
#	ax.set_ylim(-0.35,0.1)
#	ax.set_xticks([-300,-150,0,150,300])
#	ax.set_xlim(-375,375)
end
plt.savefig("PES_$(num_laps)x$(num_ions)_SOC$(tso)_tPh$(tPh)_SOCPh$(tsoPh)_multiState_cuts.pdf")
