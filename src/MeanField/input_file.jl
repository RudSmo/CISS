module InputFile

#src_path=expanduser("~/Libraries/CISS/src/MeanField/")
#include("$(src_path)hamiltonian.jl")
#include("$(src_path)iosave.jl")
#include("$(src_path)multiopenint.jl")
#include("$(src_path)classical_phonons.jl")
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
tend = 200
numstep = Int(tend/abs(dt))
writeMode = "minimal"


# Geometry Parameters
t0 = 0.04
eu = 0.65
ed = 0.65
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/40
num_laps = 1
ions = 3
l = Int(num_laps * ions)
tph = 0.1*t0
tsoph = 0.1*tso
ws = ones(l)*1.4*t0
rhoIn = "pure"
# rhoIn = "equilibrium"
phIn = "TDSCF"
# phIn = "Wigner"
popElIn = [1,1,0,0,0,0]
popPhIn = [0,0,0]
vSampl = "v1"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
qnIn,pnIn = InitialPh(popPhIn,ws,phases)

## Leads
tL = zeros(2) #[t/4,t/4] 
tR = tL #[t/4,t/4] 
V  = zeros(numstep)#[1.2 for i in 1:numstep]
nPoles = 1
beta = 1/0.025

## Central

#hcen0 = [eu 0; 0. ed] #FranssonHamiltonian(l,[eu,ed],t,tso,a,c,num_laps)


path = "data"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/L$(l)_dt$(abs(dt))"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",l),("t0",t0),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,dt,numstep,tL,tR,beta,nPoles,V,dataname2,resampling,rhoIn,phIn,vSampl,popElIn,qnIn,pnIn)
    Average(ka,dataname1,num_avg)
end

end

