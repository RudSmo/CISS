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

GammaL = zeros(2*N,2*N)*(1+0*1im)
GammaR = GammaL
GammaL[1,1] = tL[1]
GammaL[2,2] = tL[2]

GammaR[end-1,end-1] = tR[1]
GammaR[end,end] = tR[2]

Id = Matrix(I,2*N,2*N)
specFuncS = []
specFuncTu = []
specFuncTd = []
En = Array(-15:0.001:10.0)
mu = 5.0
for E in En
	GrL = inv(Id * (E+mu) - H + 1im/2 * GammaL)
	GaL = GrL'
	Lambda = GrL * GammaL * GaL
	GrR = inv(Id * (E-mu) - H + 1im/2 * GammaR)
        GaR = GrR'
	Lambda += GrR * GammaR * GaR
	push!(specFuncS,sum([Lambda[i,i]-Lambda[i+1,i+1] for i in 1:2:2*N]))
	push!(specFuncTu,sum([Lambda[i,i] for i in 1:2:2*N]))
	push!(specFuncTd,sum([Lambda[i,i] for i in 2:2:2*N]))
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

ax1.plot(En,specFuncTd,color="blue",linewidth=lw,label=L"$\sigma = \downarrow$")
ax1.plot(En,specFuncTu,color="red",linewidth=lw,label=L"$\sigma = \uparrow$")
ax1.set_ylabel(L"$\sum_i \Lambda_{i\sigma i\sigma}$")
ax1.set_xticklabels([])
ax1.tick_params(direction="in")
ax1.set_ylim(-1,200)
#ax1.set_xlim(-10,10)
ax1.legend(loc="upper left")

ax2.plot(En,specFuncTu,color="red",linewidth=lw)
ax2.plot(En,specFuncTd,color="blue",linewidth=lw)
ax2.set_ylim(0.0,0.04)
ax2.set_xlim(-1,0)
ax2.set_yticks([0,0.01,0.02,0.03,0.04])

ax3.plot(En,specFuncS,color="k",linewidth=lw)
ax3.set_ylabel(L"$\sum_i \Lambda_{i\uparrow i\uparrow} - \Lambda_{i\downarrow i\downarrow}$")
ax3.set_xlabel(L"$\varepsilon \, [t_0]$")
ax3.tick_params(direction="in")
#ax3.set_xlim(-9,0)

ax4.plot(En,specFuncS,color="k",linewidth=lw)
ax4.set_xlim(-1,0)
ax4.set_ylim(-0.001,0.001)

plt.savefig("SpectralFunction_P$(P)_V$(2*mu).pdf")



