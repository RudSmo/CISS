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
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
t0 = 0.04
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/40
num_laps = 1
ions = 1
l = Int(num_laps * ions)
tph = 0.0#0.1*t0
tsoph = 0.0#0.1*tso
ws = ones(l)*1.4*t0

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    #phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data/L3_dt0.001_1_final_Quantities"
end

popElIn = [1,0]
popPhIn = [0]
vSampl = "v1"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
qnIn,pnIn = InitialPh(popPhIn,ws,phases)

## Leads
tL = [t0/4 0;0. t0/4] #only up-up and down-down processes allowed
tR = tL #[t/4,t/4] 

spL = 0.0
spR = 0.0
V = 0.5
muLd = 0.5*(V-spL) + 6*t0
muLu = muLd+spL + 6*t0
muRd = -0.5*(V+spR) + 6*t0 
muRu = muRd+spR + 6*t0

mu = [muLu,muLd,muRu,muRd]
nPoles = 50

betaPh = 1/0.025
beta = betaPh*ones(4)



path = "data"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/L$(l)_dt$(abs(dt))_1"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",l),("t0",t0),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(l,E,t0,tso,a,c,num_laps,dt,numstep,tL,tR,beta,nPoles,mu,dataname2,rhoIn,pathIn,popElIn)
    #Average(ka,dataname1,num_avg)
end

end

