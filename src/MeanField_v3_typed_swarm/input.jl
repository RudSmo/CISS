include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")
include("structures.jl")
include("electronic_integration.jl")
using Profile
using ProfileView

using LinearAlgebra
using DelimitedFiles
# -------------- Parameters --------------------------


# Hamiltonian Parameters
num_laps = 1
num_ions = 3
N = num_laps * num_ions
ws  = 0.002*ones(N)
T0  = 0.04/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso
model = FranssonHamiltonian
modelStr = "Fransson"

# Simulation Parameters
dt = 0.0002
tend = 180 
numstep = Int(round(tend/abs(dt)))
writeMode = "minimal"


# Electronic integrators
elIntegrator = RK4_VNL!
#elIntegrator = RK45_VNL
#elIntegrator = RK45_VNL_RK4_Curs

# Classical vibrations integrators
# Non-symplectic
#clIntegrator = RK4_intertwined
#clIntegrator = RK4_sequential

# Symplectic
#clIntegrator = SymplecticEulerStep 
#clIntegrator = SymplecticVerletStep
#clIntegrator = SymplecticRuth3Step
clIntegrator = SymplecticRuth4Step


# true = propagate as usuals
# false = use initial density matrix from previous computation

initProp  = true
if initProp
#    rhoIn = "pure"
    rhoIn = "equilibrium"
    pathIn = ""
else
    rhoIn = "fromFile"
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)
popElIn[1] = 1.0

## Leads
P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)
v = 70
V = v*0.2*T0
nPoles = 15
beta = ws[1]/0.025 

xnIn = zeros(N)
pnIn = zeros(N)

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


num_avg = 100000
path = "L$(num_ions)_test/V$(v)/N$(num_avg)"
if !isdir(path)
	mkpath(path)
end


if modelStr == "Fransson"
	dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_dt$(abs(dt))"
elseif modelStr == "Holstein"
	dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
end


dataname2 = "$(dataname1)"

params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

parampath = "$(path)/paramfile.txt"
purpose = "Test of Fransson model."
SaveParams(params,parampath,purpose)

@timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,num_avg,elIntegrator,clIntegrator,popElIn)
