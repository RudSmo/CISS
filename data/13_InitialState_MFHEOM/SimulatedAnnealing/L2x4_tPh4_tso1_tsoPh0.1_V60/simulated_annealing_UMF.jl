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
xDat = readdlm("second_cycle/L2x4_tso_0.5_tPh_2.0_tsoPh_0.05_simAnn_X.dat")
xnIn,pnIn = xDat[end,:], zeros(N)
## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])
tR = tl * [1,0]#ones(2)

v = 85
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

function SimulatedAnnealing(X0,rho,model,HamPar,p_init,p_final,trial,cycle,eps)
	HamPar.X = X0
	acc = 0
	accept = false
	T_init = -1/log(p_init)
	T_final = -1/log(p_final)
	inc = (T_init-T_final)/cycle
	U_c = MFPotential(rho,HamPar)
	E_c = abs(sum(U_c))
	E_l = [E_c]
	T_l = [T_init]
	U_l = [U_c]
	X_l = [X0]
	X_c = X0
	T = T_init
	for i in 1:cycle
		for j in 1:trial
			X_i = X_c + eps*randn(HamPar.N)
			HamPar.X = X_i
			U_i = MFPotential(rho,HamPar)
			r = rand()
			E_i = real(sum(U_i))
			E_c = real(sum(U_c))
			k = exp(-(E_i-E_c)/T)
			if E_i > E_c
				if r < min(1,k)
					accept = true
				else
					accept = false
				end
			else
				accept = true
			end
			if accept == true
				E_c = E_i
				U_c = U_i
				E_c = E_i
				X_c = X_i
				HamPar.X = X_c
				acc += 1
			end
		end
		push!(E_l,E_c)
		push!(X_l,X_c)
		push!(U_l,U_c)
		T -= inc
		push!(T_l,T)			
	end
	return E_l,U_l,T_l,X_l
end

rho = read(h5open("second_cycle/L2x4_noPh_V240.0/L8_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_1_final_Quantities")["Rho"])
X0 = readdlm("second_cycle/L2x4_noPh_V240.0/L8_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_1_phonons.dat")[end,2:2:end] 
HamPar.X = X0
eps = 5*1e-3
p_init = 0.9
p_final = 0.01
trial = 200
cycle = 2000
E_l,U_l,T_l,X_l = SimulatedAnnealing(X0,rho,model,HamPar,p_init,p_final,trial,cycle,eps)
datU = zeros(length(U_l),1+length(U_l[1]))
datX = zeros(length(X_l),length(X_l[1]))
for i in 1:length(X_l)
	datU[i,1] = E_l[i]
	datU[i,2:end] = abs.(U_l[i])
	datX[i,:] = X_l[i]
end
writedlm("L$(num_laps)x$(num_ions)_tso_$(tso)_tPh_$(tPh[1])_tsoPh_$(tsoPh)_simAnn_U.dat",datU)
writedlm("L$(num_laps)x$(num_ions)_tso_$(tso)_tPh_$(tPh[1])_tsoPh_$(tsoPh)_simAnn_X.dat",datX)
writedlm("L$(num_laps)x$(num_ions)_tso_$(tso)_tPh_$(tPh[1])_tsoPh_$(tsoPh)_simAnn_T.dat",T_l)

