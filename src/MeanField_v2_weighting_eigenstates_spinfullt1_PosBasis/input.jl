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
num_laps = 3
num_ions = 3
N = num_laps * num_ions
T0  = 0.04
t0  = -T0
eu  = -0.24
ed  = -0.24
E   = [eu,ed]
Kph = 0.63  # in eV/A^2
Mph = 60014.69344356093  # in 1 / (eV A^2 )
ws  = sqrt(Kph/Mph)      # 3 meV
a   = 1.0
c   = 1.0
tso = 0.001
tPh = [0.12,-0.12]  # in eV/A
tsoPh = 0.0
model = EPeierlsHamiltonian
modelStr = "Fransson"

# Simulation Parameters
dt = 0.00005/T0
tend = 800/T0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    rhoIn = "pure"
    #rhoIn = "equilibrium"
    #phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)
popElIn[9] = 1.0
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
xnIn,pnIn = [3.06, 12.98, 0.86, -14.62, -13.78, -1.42, -2.24, 5.88],[0.368, 0.238, -0.134, -0.154, -0.054, -0.118, 0.096, 0.118]


## Leads
P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)
V  = 0.2 * T0 * 48
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
			    Kph,
			    Mph,
			    ws,
                            xnIn,
                            pnIn,
			    tL,
			    tR)



path = "L$(num_laps)x$(num_ions)_tests"
if !isdir(path)
	mkpath(path)
end


num_avg = 1
for ka in 1:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
    end
     
     
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
#    Average(ka,dataname1,num_avg)
end

end

