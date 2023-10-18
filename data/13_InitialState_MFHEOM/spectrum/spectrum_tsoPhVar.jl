src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

src=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_swarm_PVar/")
include("$(src)classical_phonons.jl")

using PyPlot
using LinearAlgebra
using DelimitedFiles
using HDF5

# Hamiltonian Parameters
num_laps = 4
num_ions = 4
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
tso = T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"
# Simulation Parameters
dt = 0.0002
tend = 3
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"


popElIn = zeros(2*N)
popElIn[1] = 0.0
popPhIn = zeros(N)

xDat = readdlm("phonon.dat")
#xnIn,pnIn = xDat[end,2:2:end],zeros(N)
#xDat[end,:], zeros(N)
num_avg = 1
beta = ws[1]/0.025
xnIn,pnIn = SampleWigner_Multi(beta,ws,N,num_avg)
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 75
V = v * 0.2 * T0
nPoles = 15

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

tso_rng = 0.0:0.1*tso:50*tso
EW = []
for tsos in tso_rng
	HamPar.SOCPh = tsos
	H = model(HamPar)
	ew = eigvals(H)
	push!(EW,ew)
end
plt.plot(tso_rng .* ws[1] * 1000 ,EW .* ws[1] * 1000)
plt.xlabel(L"$\lambda_1 \, [\mathrm{meV}]$")
plt.ylabel(L"$E\,[\mathrm{meV}]$")
plt.savefig("L$(num_laps)x$(num_ions)_tsoPh_$(tsoPh)_X_zeros_tsoPhVar.pdf")

