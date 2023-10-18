src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

using PyPlot
using LinearAlgebra
using DelimitedFiles
using HDF5

# Hamiltonian Parameters
num_laps = 2
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
xnIn,pnIn = xDat[1,2:2:end],zeros(N)
#xnIn,pnIn = zeros(N),zeros(N)
#xDat[end,:], zeros(N)
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 75
V = v * 0.2 * T0
nPoles = 15
beta = ws[1]/0.025

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

ES_pop = []
fac_t = 1/(0.002 * 1e3 * 4.135668)
time = xDat[2:end,1]
state = N 
for i in 2:length(time)+1
	HamPar.X = xDat[i,2:2:end]
	H = model(HamPar)
	ev = eigvecs(H)[state,:]
	pop = [real(ev[j:j+1]' * [1 0; 0 1.] * ev[j:j+1]) for j in 1:2:2*N]	
	push!(ES_pop,pop)
end
plt.plot(time * fac_t,ES_pop)
plt.xlabel(L"$t \, [\mathrm{ps}]$")
plt.ylabel(L"$n_j$")
plt.savefig("L$(num_laps)x$(num_ions)_tso_$(tso)_tsoPh_$(tsoPh)_XVar_pop_es_$(state).pdf")

