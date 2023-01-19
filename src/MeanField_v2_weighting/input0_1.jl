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
num_laps = 4
num_ions = 4
N = num_laps * num_ions
t0  = 0.04
eu  = -6*t0
ed  = -6*t0
E   = [eu,ed]
ws  = 0.05*t0*ones(N)
gPh = [0.0,0.0]
a   = 1.0
c   = 1.0
tso = t0/40
tPh = 0.0 #0.1*t0
tsoPh = 0.0 #0.1*tso
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
#model = FranssonHamiltonianNH
modelStr = "Fransson"


# Simulation Parameters
dt = 0.002/t0
tend = 250/t0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    #rhoIn = "pure"
    rhoIn = "equilibrium"
    phIn = "TDSCF"
    #phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data_single/L1_Delta$(eu-ed)_gph$(gph)_ws$(ws)_dt$(abs(dt))_final_Quantities"
end
popElIn = 1.0*[1,1,0,0,0,0,0,0]
popPhIn = [0,0,0,0,0,0,0,0]
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
x0 = 1.0
xnIn,pnIn = [-4.28, 13.0, -2.66, -6.68, 4.82, 12.32, -0.66, 9.58, -3.98, -4.06, 12.22, 12.34, 7.76, 3.52, 5.96, -11.82], 
 	    [-0.022, -0.082, -0.024, 0.132, -0.128, 0.188, -0.014, 0.168, 0.016, 0.02, -0.002, -0.048, -0.204, -0.034, -0.012, 0.034] 


## Leads
P = 0.5
tL = t0/4 * (ones(2) .+ P * [1,-1])  
tR = t0/4 * ones(2)
V  = 0*t0
nPoles = 20
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



path = "Tests_4x4/V0.48/KM"#WignerTraj/Fransson"
if !isdir(path)
	mkpath(path)
end


num_avg = 2
for ka in 2:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh$(tPh)_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_X0_-10_dt$(abs(dt))"
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

