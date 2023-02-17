module InputFile

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_spins.jl")
include("structures.jl")
using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using .ClassicalSpins
using .CustomStructures
# -------------- Parameters --------------------------


# Hamiltonian Parameters
N   = 4
t0  = 1.0
Jsd = 0.5
eu  = 0.0#-6*t0
ed  = 0.0#-6*t0
E   = [eu,ed]
Hext = [[0,0,0.],[0,0,0.],[0,0,0.]]
#model = Hamiltonian_Kondo
#modelStr = "Kondo"

model = Hamiltonian_KondoPeierls
modelStr = "Kondo-Peierls"

# Simulation Parameters
dt = 0.002/t0
tend = 2000/t0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
popElIn = zeros(2*N)
function RandomInitialMag()
	v = randn(3)
	return v/norm(v)	
end
pos = [1,2,3]
magIn   = [RandomInitialMag() for i in 1:length(pos)]

## Leads
P = 0.5
tL = (ones(2) .+ P * [1,-1])  
tR = ones(2)
V  = 1.5
nPoles = 30
beta = 1/0.025
HamPar = HamiltonianParams( N,
                            model,
                            modelStr,
                            E,
                            t0,
			    Jsd,
                            magIn,
			    pos,
			    Hext,
			    tL,
			    tR)


if modelStr == "Kondo"
	path = "Kondo_test"
elseif modelStr == "Kondo-Peierls"
	path = "KondoPeierls_test"
end
if !isdir(path)
	mkpath(path)
end

dataname = "$(path)/L$(N)_Jsd$(Jsd)_V$(V)_P$(P)_multi_dt$(abs(dt))"
params = [("dt",dt),("Nstep",numstep),("N",1),("tL",tL),("tR",tR),("Jsd",Jsd),("V",V),("nPoles",nPoles),("beta",beta)]
parampath = "$(path)/paramfile.txt"
purpose = "Test of Spin Peierls model."
IOsave.SaveParams(params,parampath,purpose)
@timev Propagate(HamPar,dt,numstep,beta,nPoles,V,dataname,rhoIn,magIn,popElIn)

end

