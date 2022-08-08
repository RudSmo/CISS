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
tend = 1000
numstep = Int(tend/abs(dt))
writeMode = "minimal"


# Geometry Parameters
eu  = 1.0
ed  = -eu
E   = [eu,ed]
ws  = 1.4*eu
gph = 0 * 0.1

#rhoIn = "pure"
rhoIn = "equilibrium"

# phIn = "TDSCF"
phIn = "Wigner"

#popElIn = [1,1,0,0,0,0]
#popPhIn = [0,0,0]
vSampl = "v1"
#phases = [rand()*2*pi for i in 1:length(popPhIn)]
#qnIn,pnIn = InitialPh(popPhIn,ws,phases)

## Leads
tL = [0.0,0.0] 
tR = [eu,eu] 
V  = 0 * 0.1 * ones(numstep)
nPoles = 15
beta = 1/0.025

path = "data_gph$(gph)_ws$(ws)"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/L1_gph$(gph)_ws$(ws)_dt$(abs(dt))"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gph",gph),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev PropagateSingle(E,ws,gph,dt,numstep,tL,tR,beta,nPoles,V,dataname2,resampling,rhoIn,phIn,vSampl)
end

end

