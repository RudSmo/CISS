module InputFile

src_path=expanduser("~/Libraries/CISS/src/MeanField_v2_weighting_eigenstates_spinfullt1/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)iosave.jl")
include("$(src_path)multiopenint.jl")
include("$(src_path)classical_phonons.jl")
include("$(src_path)structures.jl")
using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
using .ClassicalPhonons
using .CustomStructures
# -------------- Parameters --------------------------

# Simulation Parameters
T0 = 0.04
t0 = -T0
dt = 0.002/T0
tend = 5000/T0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Hamiltonian Parameters
num_laps = 1
num_ions = 50
N = num_laps * num_ions
eu  = 0.0
ed  = 0.0
E   = [eu,ed]
ws  = 0.002*ones(N)
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.0
tPh = [-0.004,-0.004]
tsoPh = 0.0
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"

if initProp
    #rhoIn = "pure"
    rhoIn = "equilibrium"
    #phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
popElIn = zeros(2*N)
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
xnIn,pnIn = InitialPh(popPhIn,ws,phases)


## Leads
P = 0.0

tL = T0/4 * (ones(2) .+ P * [1,-1])
tR = T0/4 * ones(2)
V = 4*T0/3
nPoles = 10
beta = 1/0.025
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
                            pnIn,
			    tL,
			    tR)



path = "data0"
if !isdir(path)
	mkdir(path)
end

num_avg = 1
for ka in 1:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPhU_$(tPh[1])_tPhD_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_dt$(abs(dt))"
    end
     
     
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
    Average(ka,dataname1,num_avg)
end

end

