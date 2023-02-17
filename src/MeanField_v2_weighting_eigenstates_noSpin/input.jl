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
num_ions = 2
N = num_laps * num_ions
t0  = 2.0
E   = 0.0
ws  = 0.1*ones(N)
gPh = 0.0*t0
tPh = 0.3
model = HolsteinHopping
modelStr = "Holstein"

# Simulation Parameters
dt = 0.001/t0
tend = 2000/t0
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
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
x0 = 0.0
xnIn,pnIn = [0.25,0.0],[-0.05,0.0]


## Leads
tL = 0.01
tR = 0.01
nPoles = 10
beta = 1/0.005
V = 3.0
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


path = "Peierls_L2_SimpleParams"
if !isdir(path)
	mkpath(path)
end

num_avg = 1
for ka in 1:num_avg
    dataname1 = "$(path)/L$(N)_gPh$(gPh)_tPh$(tPh)_V$(V)_w0$(ws[1])_dt$(abs(dt))" 
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Holstein model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
    #Average(ka,dataname1,num_avg)
end

end

