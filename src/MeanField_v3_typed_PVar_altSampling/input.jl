include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")
include("structures.jl")
include("electronic_integration.jl")

using LinearAlgebra
using DelimitedFiles
# -------------- Parameters --------------------------


# Hamiltonian Parameters
num_laps = 1
num_ions = 3
N = num_laps * num_ions
ws  = 0.002*ones(N)
Hop = 0.04
T0  = Hop/ws[1]  # in eV
t0  = -T0
eu  = -6*T0
ed  = -6*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.0 * T0/40
tPh = 1.0*T0*ones(2)
tsoPh = 0.0*tso
model = FranssonHamiltonian
modelStr = "Fransson"


# Simulation Parameters
ntraj = 10
dt = 0.0002
tend = 4
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


rhoIn = "pure"
#    rhoIn = "equilibrium"
#    phIn = "TDSCF"
phIn = "Wigner"


popElIn = zeros(2*N)
popElIn[1] = 0.0
popPhIn = zeros(N)
xnIn,pnIn = zeros(N,ntraj),zeros(N,ntraj)

## Leads
P = 0.0
tl = T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])  
tR = tl * [1,0]#ones(2)

v = 85
V = v * 0.2 * T0
nPoles = 15
beta = ws[1]/0.025 

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


version = 0
path = "L$(num_laps)x$(num_ions)_Peierls_V$(V)"
if !isdir(path)
	mkpath(path)
end


if modelStr == "Fransson"
	dataname = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_dt$(abs(dt))"
elseif modelStr == "Holstein"
	dataname = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
end
     
     

params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

parampath = "$(path)/paramfile.txt"
purpose = "Test of Fransson model."
SaveParams(params,parampath,purpose)
    
if phIn == "Wigner"
	X0,P0 = SampleWigner_Multi_ntraj(beta,ws,N,ntraj)
	HamPar.X[:,:]  = X0 
	HamPar.P[:,:]  = P0
end
@timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname,rhoIn,elIntegrator,clIntegrator,popElIn)

