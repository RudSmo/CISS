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
tsoPh = 5*tso
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

xDat = readdlm("second_cycle/L2x4_tso_0.5_tPh_2.0_tsoPh_0.05_simAnn_X.dat")
xnIn,pnIn = xDat[end,:], zeros(N)
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

function PhononHamiltonian(HamPar)
        H = zeros(2*HamPar.N,2*HamPar.N)
        for i in 1:HamPar.N
                H[2*i-1,2*i-1] = 0.5* HamPar.X[i]^2
                H[2*i,2*i] = 0.5*HamPar.X[i]^2
        end
        return H
end

function MFPotential(rho,HamPar)
	UMF = []
	for mode in 1:length(HamPar.X)
		Uvib = 0.5 * HamPar.X[mode]^2 
		Helvib = FranssonHamiltonian_Mode(HamPar,mode) 
		Uelvib = tr(Helvib * rho)
		push!(UMF,Uelvib+Uvib)
	end
	return UMF	
end

rho = read(h5open("initial_prop_V85/L2x4_noPh_V340.0/L8_tPh_0.0_thPd_0.0_SOC0.5_SOCPh0.0_V340.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_1_final_Quantities")["Rho"])
X_rng = -5:0.1:5

for site in 1:N
	UMF = []
	for X in X_rng
		HamPar.X[site] = X
		umf = MFPotential(rho,HamPar)[site]
		push!(UMF,real(umf))
	end
	println(X_rng[findmin(UMF)[2]])
	plt.plot(X_rng,UMF,label=site)
end
plt.legend()
plt.show()
