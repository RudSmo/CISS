module InputFile

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("structures.jl")
using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using .CustomStructures
# -------------- Parameters --------------------------


# Hamiltonian Parameters
num_laps = 1 
num_ions = 50
N = num_laps * num_ions
T0  = 0.04
t0  = -T0
eu  = 0.0 #-6*T0
ed  = 0.0 #-6*T0
E   = [eu,ed]
a   = 1.0
c   = 1.0
tso = T0/40
modelStr = "Fransson"
model = FranssonHamiltonian

# Simulation Parameters
dt = 0.002/T0
tend = 500/T0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    #rhoIn = "pure"
    rhoIn = "equilibrium"
    pathIn = ""
else
    rhoIn = "fromFile"
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)


## Leads
P = 0.0
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)

Vtype = "const"
Vdrop = "sym"
path = "L$(N)_V$(Vtype)_var"
if !isdir(path)
	mkpath(path)
end

for vs in Array(0.0:T0:3*T0)

	path = "L$(N)_V$(Vtype)_var/V$(vs)"
	if !isdir(path)
		mkpath(path)
	end
	
	V = (vs,)
	#V = (80/T0,7.5/T0) # (t0,tau)
	#Vtype = "lorentz" 
	#Vdrop = "left"

	nPoles = 10
	beta = 1/0.025 
	HamPar = HamiltonianParams( num_laps,
				    num_ions,
				    N,
				    model,
				    modelStr,
				    E,
				    a,
				    c,
				    t0,
				    tso,
				    tL,
				    tR)






	if modelStr == "Fransson"
		if Vtype == "const"
			dataname = "$(path)/L$(N)_SOC$(tso)_V$(V[1])_P$(P)_noEps_dt$(abs(dt))"
		elseif Vtype == "lorentz"
			dataname = "$(path)/L$(N)_SOC$(tso)_Vt0$(V[1])_Vtau$(V[2])_P$(P)_noEps_dt$(abs(dt))"
		end
	end



	params = [("dt",dt),("Nstep",numstep),("N",1),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

	parampath = "$(path)/paramfile.txt"
	purpose = "Test of Fransson model."
	IOsave.SaveParams(params,parampath,purpose)

	Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,Vtype,Vdrop,dataname,rhoIn,popElIn)
end

end

