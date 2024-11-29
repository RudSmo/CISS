src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
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
num_laps = 1
num_ions = 5
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
P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
Gamma_L,Gamma_R = GammaMatrices1d(N,tL,tR)
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

P = "P"
filL = open("gammaL_$(P)_L$(num_laps)x$(num_ions).dat","w")
filR = open("gammaR_$(P)_L$(num_laps)x$(num_ions).dat","w")
datPh = readdlm("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat")
X = datPh[2:end,2:2:end]
t = datPh[2:end,1]
for i in 1:length(X[:,1])
	hp = HamiltonianParams( num_laps,
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
	H = model(hp)
	ew,U = eigen(H)
	Gl = U' * Gamma_L * U
	Gr = U' * Gamma_R * U
	write(filL,"$(t[i]) ")
	write(filR,"$(t[i]) ")
	for j in 1:length(Gl[1,:])
		write(filL,"$(real(Gl[j,j])) ")
	end
	for j in 1:length(Gr[1,:])
                write(filR,"$(real(Gr[j,j])) ")
        end	
	write(filL,"\n")
	write(filR,"\n")
end
close(filL)
close(filR)
