using LinearAlgebra
using PyPlot
src_path=expanduser("/Users/rudolfsmorka/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

using .Hamiltonian
using .CustomStructures

num_laps = 2
num_ions = 4
N = num_laps * num_ions
t0  = 1.0
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.0
tsoPh = 0.0
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn = zeros(N)
pnIn = zeros(N)

P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])
tR = t0/4 * ones(2)
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
                            xnIn,
                            pnIn,
                            tL,
                            tR)


H = model(HamPar)

Gamma = zeros(2*N,2*N)*(1+0*1im)
Gamma[1,1] = tL[1]
Gamma[2,2] = tL[2]
#Gamma[end-1,end-1] = tR[1]
#Gamma[end,end] = tR[2]

Id = Matrix(I,2*N,2*N)
specFuncS = []
specFuncTu = []
specFuncTd = []
En = Array(-10:0.001:1.0)
site = 7
s = 2*site - 1

for E in En
	Gr = inv(Id * E - H + 1im/2 * Gamma)
	Ga = Gr'
	Lambda = Gr * Gamma * Ga
	push!(specFuncS,Lambda[s,s]-Lambda[s+1,s+1])
	push!(specFuncTu,Lambda[s,s])
	push!(specFuncTd,Lambda[s+1,s+1])
end

fig = plt.figure(figsize=(10,7))
px = 0.09
py = 0.55
sx = 0.87
sy = 0.42
dy = 0.03
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([0.63,0.63,0.3,0.27])
ax3 = plt.axes([px,py-dy-sy,sx,sy])
ax4 = plt.axes([0.63,0.15,0.3,0.27])

lw = 3
ax1.plot(En,specFuncTu,color="red",linewidth=lw,label=L"$\sigma = \uparrow$")
ax1.plot(En,specFuncTd,color="blue",linewidth=lw,label=L"$\sigma = \downarrow$")
ax1.set_ylabel(L"$\Lambda_{1\sigma 1\sigma}$")
ax1.set_xticklabels([])
ax1.tick_params(direction="in")
ax1.set_ylim(-1,40)
ax1.set_xlim(-9,0)
ax1.legend(loc="upper left")

ax2.plot(En,specFuncTu,color="red",linewidth=lw)
ax2.plot(En,specFuncTd,color="blue",linewidth=lw)
#ax2.set_ylim(0.0,0.02)
ax2.set_ylim(-0.01,0.01)
ax2.set_xlim(-1,0)

ax3.plot(En,specFuncS,color="k",linewidth=lw)
ax3.set_ylabel(L"$\Lambda_{1\uparrow 1\uparrow} - \Lambda_{1\downarrow 1\downarrow}$")
ax3.set_xlabel(L"$\varepsilon \, [t_0]$")
ax3.tick_params(direction="in")
ax3.set_xlim(-9,0)

ax4.plot(En,specFuncS,color="k",linewidth=lw)
ax4.set_xlim(-1,0)
#ax4.set_ylim(0.0,0.02)
ax4.set_ylim(-0.01,0.01)

plt.savefig("SpectralFunction_P0.5_site$(site).pdf")



