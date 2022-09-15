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
tend = 120
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
num_laps = 1
ions = 40
l = Int(num_laps * ions)
t0 = ones(2*l-2)
eu = 1.0
ed = -1.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = zeros(l)
tso[Int(l/2)] = 1.0
tso[Int(l/2)+1] = 1.0
tso[Int(l/2)+2] = 1.0
tso[Int(l/2)+3] = 1.0
tso[Int(l/2)-1] = 1.0
tso[Int(l/2)-2] = 1.0
tso[Int(l/2)-3] = 1.0
tph = zeros(l)   
tsoph = zeros(l)
ws = 1.4*t0[1]*ones(l)

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    phIn = "TDSCF"
    #phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data/L3_dt0.001_1_final_Quantities"
end

popElIn = zeros(2*l)
popElIn[1] = 1.0
popElIn[2] = 1.0

popPhIn = zeros(l)
vSampl = "v1"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
qnIn,pnIn = InitialPh(popPhIn,ws,phases)

## Leads
tL = 0.0 * [1.0,1.0] 
tR = tL #[t/4,t/4] 
V  = 0.0*ones(numstep)
nPoles = 50
beta = 1/0.025


path = "data_comparison"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/L$(l)_Hop+Onsite+SOC_dt$(abs(dt))"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",l),("t0",t0),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,dt,numstep,tL,tR,beta,nPoles,V,dataname2,resampling,rhoIn,phIn,vSampl,pathIn,popElIn,qnIn,pnIn)
    #Average(ka,dataname1,num_avg)
end

end

