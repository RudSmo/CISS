src_path=expanduser("~/Libraries/CISS/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)iosave.jl")
include("$(src_path)multiopenint.jl")
include("$(src_path)classical_phonons.jl")
include("$(src_path)structures.jl")
include("$(src_path)electronic_integration.jl")

using LinearAlgebra
using DelimitedFiles
using PyPlot
# Hamiltonian Parameters
num_laps = 7
num_ions = 3
N = num_laps * num_ions
ws  = 0.2*0.002*ones(N)
T0  = 0.04/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = 0.1 * T0 * ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"

# Simulation Parameters
dt = 0.0003
tend = 100
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"

xnIn,pnIn = zeros(N), zeros(N)
## Leads
P = -0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
V = 60 * 0.2 * T0
nPoles = 15
beta = ws[1]/0.025
initProp  = true

popElIn = zeros(2*N)
popPhIn = zeros(N)
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

datPh = readdlm("../0/phonons_avg_n1001_M.dat")
X = datPh[2:40:end,2:2:end]
ew = []
for i in 1:length(X[:,1])
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
                            X[i,:],
                            pnIn,
                            tL,
                            tR)
        H = model(HamPar)
        push!(ew,eigvals(H))
end

plt.plot(datPh[2:40:end,1],ew)
plt.show()
