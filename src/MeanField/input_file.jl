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
dt = 0.002
tend = 25
numstep = Int(tend/abs(dt))
writeMode = "minimal"


# Geometry Parameters
t = 0.04
eu = 0.65
ed = 0.65
E = [eu,ed]
a = 1.0
c = 1.0
tso = t/40
num_laps = 1
ions = 2
l = Int(num_laps * ions)
tph = 0.1*t
tsoph = 0.1*tso
ws = ones(l)*1.4*t


## Leads
tL = [t/4,t/4] 
tR = [t/4,t/4] 
V  = [1.2 for i in 1:numstep]
nPoles = 40
beta = 1/0.025

## Central

#hcen0 = [eu 0; 0. ed] #FranssonHamiltonian(l,[eu,ed],t,tso,a,c,num_laps)


path = "test"
mkpath(path)
num_avg = 5
resampling = true
rsmpl_stp = Int(5/dt)
for ka in 1:num_avg
    dataname1 = "$(path)/L$(l)_V$(V[1])_dt$(abs(dt))"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",l),("t",t),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(l,E,t,tso,a,c,num_laps,tph,tsoph,ws,dt,numstep,tL,tR,beta,nPoles,V,dataname2,resampling,rsmpl_stp)
    Average(ka,dataname1,num_avg)

end

end

