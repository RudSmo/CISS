module InputFile

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")
include("structures.jl")
using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using .ClassicalPhonons
using .CustomStructures
# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 400
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Hamiltonian Parameters
num_laps = 1
num_ions = 2
N = num_laps * num_ions
t0  = 1.0
eu  = 0.0
ed  = 0.0
E   = [eu,ed]
ws  = 0.1*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = 0.0
tPh = 0.0
tsoPh = 0.0
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    phIn = "TDSCF"
    #phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
popElIn = [1.0,0,0.0,0]
popPhIn = [0,0]
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
xnIn,pnIn = [1.0,0],[0.0,0]#InitialPh(popPhIn,ws,phases)


HamPar = HamiltonianParams( num_laps,
                            num_ions,
                            N,
                            model,
                            modelStr,
                            E,
                            a,
                            c,
                            t0,
                            tPh,
                            tso,
                            tsoPh,
                            gPh,
                            ws,
                            xnIn,
                            pnIn)


## Leads
P = 0.0
tL = 0.0 * (ones(2) .+ P * [1,-1])  
tR = 0.0 * ones(2)
V  = 0.0
nPoles = 1
beta = 1/0.025

path = "TwoSite_Comp"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
for ka in 1:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_w0$(ws[1])_P$(P)_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_dt$(abs(dt))"
    end
     
     
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
    #Average(ka,dataname1,num_avg)
end

end

