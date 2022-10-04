module InputFile

#src_path=expanduser("~/Libraries/CISS/src/MeanField/")
#include("$(src_path)hamiltonian.jl")
#include("$(src_path)iosave.jl")
#include("$(src_path)multiopenint.jl")
#include("$(src_path)classical_phonons.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")

using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 500.0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
E = [zeros(2),zeros(2),zeros(2)];
t0 = [1.0,1.0,1.0];
gSL = [[0.25 0; 0 0.25],[0.25 0; 0 0.25]];
tso = [0.0,t0[2]/4,0.0]
num_laps = [1,3,1]
ions = [15,4,15]
l = Int.([n*i for (n,i) in zip(num_laps,ions)]);
a = ones(3);
c = a
tph = zeros(3)
tsoph = zeros(3)
ws = [zeros(i) for i in l]
qs = ws
H = Hamiltonian.FranssonHamiltonianExtended(l,E,t0,gSL,tso,a,c,num_laps,tph,tsoph,ws,qs)
if initProp
    #rhoIn = "pure"
    rhoIn = "equilibriumClosed"
    phIn = "TDSCF"
    # phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data/L3_dt0.001_1_final_Quantities"
end

popElIn = [1,1,0,0,0,0,0,0,0,0]

## Leads
tL = 0.25*ones(2) 
tR = tL #[t/4,t/4] 
V  = 0.0  #[1.2 for i in 1:numstep]
muL = V/2 #+ 6*t0
muR = -V/2 #+ 6*t0
nPoles = 2
beta = 1/0.025


path = "data/Open_Extended"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
resampling = false
rsmpl_stp = 15
for ka in 1:num_avg
    dataname1 = "$(path)/Le$(l[1])_LS$(l[2])_laps$(num_laps[2])_ions$(ions[2])_V$(V)_tso$(tso[2])_gSL$(gSL[1][1])_dt$(abs(dt))_2"
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",l),("t0",t0),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(Int(sum(l)),H,dt,numstep,tL,tR,beta,nPoles,muL,muR,dataname2,rhoIn,popElIn)
    #Average(ka,dataname1,num_avg)
end

end

