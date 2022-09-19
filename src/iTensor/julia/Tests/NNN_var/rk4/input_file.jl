module InputFile

src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_el/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)iosave.jl")
include("$(src_path)multiopenint.jl")

using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 50
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
t1 = 1.0
t2 = 1.0
tns = [t1,t2]
l = 4
H = Hamiltonian.NNN(l,tns)
if initProp
    rhoIn = "pure"
    # rhoIn = "equilibrium"
    phIn = "TDSCF"
    # phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data/L3_dt0.001_1_final_Quantities"
end
popElIn = zeros(2*l)
for i in 1:2
	popElIn[i] = 1
end
#popElIn = [1,1,0,0,0,0,0,0]

## Leads
tL = zeros(2) 
tR = tL #[t/4,t/4] 
V  = 0.0  #[1.2 for i in 1:numstep]
muL = V/2 
muR = -V/2 
nPoles = 1
beta = 1/0.025

dataname = "L$(l)_n$(length(tns))_dt$(abs(dt))"

params = [("dt",dt),("Nstep",numstep),("N",l),("tns",tns),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

parampath = "paramfile.txt"
purpose = "Test of Fransson model."
IOsave.SaveParams(params,parampath,purpose)

@timev Propagate(l,H,dt,numstep,tL,tR,beta,nPoles,muL,muR,dataname,rhoIn,popElIn)

end

