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
N = 21
dt = 0.001
tend = 100.0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = false
# true = propagate as usuals
# false = use initial density matrix from previous computation
rhoIn = "pure"
popIn = zeros(N)
popIn[Int(ceil(N/2))] = 1



# Hoppings
t0 = 1.0
h = 1000.0
tr = t0
tl = t0/h 
H = Hamiltonian.HatanoNelson(N,tr,tl) 

## Leads
tL = 0.0 
tR = tL 
V  = 1.0
nPoles = 1
beta = 1/0.025


path = "N$(N)_closed_h$(h)"
if !isdir(path)
	mkdir(path)
end

dataname = "$(path)/L$(N)"
params = [("dt",dt),("Nstep",numstep),("N",N),("tl",tl),("tr",tr),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]
parampath = "$(path)/paramfile.txt"
purpose = "Test of Fransson model."
IOsave.SaveParams(params,parampath,purpose)
Propagate(N,H,dt,numstep,tL,tR,beta,nPoles,V,dataname,rhoIn,popIn)

end

