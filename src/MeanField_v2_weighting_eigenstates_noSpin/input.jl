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


# Hamiltonian Parameters
num_laps = 1
num_ions = 3
N = num_laps * num_ions
T0  = 0.04
t0  = -T0
E   = -0.24
ws  = 0.01*ones(N)
gPh = 0.0*T0
tPh = 0.0
model = HolsteinHopping
modelStr = "Holstein"

# Simulation Parameters
dt = 0.002/T0
tend = 600/T0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    phIn = "TDSCF"
#    phIn = "Wigner"
#    phIn = "unocc"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
popElIn = zeros(N)
popElIn[1] = 1.0
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
x0 = 0.0
xnIn,pnIn = [1,-1.,0.5],[0.,0,0]
#zeros(N),zeros(N)
#[-0.1,0.3,0.05,0],[0.002,-0.01,0.0005,0]


## Leads
tL = 0.004
tR = 0.004
nPoles = 10
beta = 1/0.025
V = 0.4
path = "L3_Leads"
if !isdir(path)
	mkpath(path)
end
HamPar = HamiltonianParams( num_laps,
       	                    num_ions,
                            N,
                            model,
                            modelStr,
                            E,
                            t0,
                            tPh,
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
			    tL,
			    tR)


num_avg = 1
for ka in 1:num_avg
    dataname1 = "$(path)/L$(N)_gPh$(gPh)_tPh$(tPh)_V$(V)_w0$(ws[1])_E$(E)_dt$(abs(dt))" 
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Holstein model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
    #Average(ka,dataname1,num_avg)
end
end

