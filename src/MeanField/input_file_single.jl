module InputFile

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")

using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using .ClassicalPhonons
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
eu  = 0.2
ed  = -eu
E   = [eu,ed]
ws  = 1.0
gph = 0.5

if initProp
    #rhoIn = "pure"
    rhoIn = "equilibrium"
    # phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
#popElIn = [1,1,0,0,0,0]
#popPhIn = [0,0,0]
vSampl = "v1"
#phases = [rand()*2*pi for i in 1:length(popPhIn)]
#qnIn,pnIn = InitialPh(popPhIn,ws,phases)

## Leads
tL = [0.0,0.0] 
tR = ones(2)
V  = 0 * 0.1 
nPoles = 15
beta = 1/0.025

path = "data_single"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gph",gph),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev PropagateSingle(E,ws,gph,dt,numstep,tL,tR,beta,nPoles,V,dataname2,resampling,rhoIn,phIn,vSampl,pathIn)
end

end

