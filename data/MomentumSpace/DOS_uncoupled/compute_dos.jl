path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting_eigenstates_spinfullt1/")
include("$(path)hamiltonian.jl")
include("$(path)structures.jl")
include("$(path)classical_phonons.jl")

using LinearAlgebra
using .Hamiltonian
using .CustomStructures
using .ClassicalPhonons

# Hamiltonian Parameters
num_laps = 1
num_ions = 50
N = num_laps * num_ions
T0  = 0.04
t0  = -T0
eu  = 0.0
ed  = 0.0
E   = [eu,ed]
ws  = 0.002*ones(N)
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.0
tPh = [-0.004,0.004]
tsoPh = 0.0
model = FranssonHamiltonian
modelStr = "Fransson"

## Leads
P = 0.0
tl = 0.004
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
beta = 1/0.025

numSample = 1000
fil = open("DOS_N_$(N)_tso_$(tso)_tPhU_$(tPh[1])_tPhD_$(tPh[2])_tsoPh_$(tsoPh)_beta_$(beta)_nS_$(numSample)_uncoupled.dat","w")

for i in 1:numSample
	xnIn,pnIn = SampleWigner_Multi(beta,ws,N) 
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
	ew = eigvals(H)
	write(fil,"$(i) ")
	for j in ew
		write(fil,"$(j) ")
	end
	write(fil,"\n")
end
close(fil)
	
