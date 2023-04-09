path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v2_weighting_eigenstates_spinfullt1/")
include("$(path)hamiltonian.jl")
include("$(path)structures.jl")
include("$(path)classical_phonons.jl")

using LinearAlgebra
using DelimitedFiles
using .Hamiltonian
using .CustomStructures
using .ClassicalPhonons


# Hamiltonian Parameters
num_laps = 3
num_ions = 4
N = num_laps * num_ions
T0  = 0.04
t0  = -T0
eu  = -0.24
ed  = -0.24
E   = [eu,ed]
ws  = 0.002*ones(N)
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = [-0.004,0.004]
tsoPh = 0.1 * tso
model = FranssonHamiltonian
modelStr = "Fransson"

## Leads
P = 0.0
tl = 0.004
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
beta = 1/0.025

Erng = -9*T0:0.00001:-3*T0

outfil = open("DOS_N_$(N)_tso_$(tso)_tPhU_$(tPh[1])_tPhD_$(tPh[2])_tsoPh_$(tsoPh)_beta_$(beta)_dynamics_uncoupled.dat","w")
phfil  = readdlm("phonons_avg_n20.dat") 
ntime  = length(phfil[2:end,1])

function DOS(N,eigv,Erng,sigma = 1e-7)
        return [1/(sqrt(2*pi)*sigma*N)*sum([exp(-(E-i)^2/(2*sigma^2)) for i in eigv]) for E in Erng]
end



for i in 2:ntime
	xnIn = phfil[i,2:2:end]
	pnIn = phfil[i,2:2:end]
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
	#dos = DOS(N,ew,Erng)
	write(outfil,"$(phfil[i,1]) ")
	for j in ew
		write(outfil,"$(j) ")
	end
	write(outfil,"\n")
end
close(outfil)
