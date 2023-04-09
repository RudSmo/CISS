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
num_laps = 6
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
tPh = [0.0,0.0]
tsoPh = 0.0
model = FranssonHamiltonian
modelStr = "Fransson"

## Leads
P = 0.0
tl = 0.004
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
GL = zeros(2*N,2*N)
GL[1,1] = tL[1]
GL[2,2] = tL[2]
GR = zeros(2*N,2*N)
GR[end-1,end-1] = tR[1]
GR[end,end] = tR[2]
beta = 1/0.025

numSample = 1
filname = "GF_DOS_N_$(N)_tso_$(tso)_tPhU_$(tPh[1])_tPhD_$(tPh[2])_tsoPh_$(tsoPh)_beta_$(beta)_tl_$(tl)_nS_$(numSample)_P_$(P)_coupled.dat"


Erng = -9*T0:0.00001:-3*T0
Sigma = -1im/2 * (GL + GR)
LDOS_ens = zeros(length(Erng),2)
LDOS_ens[:,1] = Erng
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
	id = Matrix(I,size(H))
	LDOS = zeros(length(Erng))
	for (i,E) in enumerate(Erng)
		Gr = inv(E*id - H + Sigma)
		As = 1im * (Gr - Gr')
		ldos = -1/(2*pi) * sum([As[i,i] for i in 1:size(As)[1]])
		LDOS[i] = ldos
	end	
	LDOS_ens[:,2] .+= LDOS

	fil = open(filname,"w")
	writedlm(fil,LDOS_ens)
	close(fil)
end
	
