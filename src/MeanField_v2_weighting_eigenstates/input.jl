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
t0  = 0.04
eu  = 0.0#-6*t0
ed  = 0.0 #-6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh = ones(2) * 0.1 * t0
a   = 1.0
c   = 1.0
tso = 0.0 # t0 / 40
tPh = 0.0 * t0 
tsoPh = 0.0*tso
model = HolsteinHopping
modelStr = "Holstein"
#model = FranssonHamiltonian
#modelStr = "Fransson"
#model = FranssonHamiltonianNH

# Simulation Parameters
dt = 0.002/t0
tend = 750/t0
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
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
x0 = 0.0
xnIn,pnIn = zeros(N),zeros(N)


## Leads
P = 0.0
tl = t0/10
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)
V  = 0*t0
nPoles = 20
beta = 1/0.0125 #1/0.025
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



#path = "Peierls_L2/Eps"
#path = "Peierls_L2/noEps_largerGap"
#path = "Holstein_L2/noEps"
path = "Holstein_L2/noEps_largerGap"
if !isdir(path)
	mkpath(path)
end


num_avg = 1
for ka in 1:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_noEps_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
    end
     
     
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    IOsave.SaveParams(params,parampath,purpose)

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
#    Average(ka,dataname1,num_avg)
end

end

