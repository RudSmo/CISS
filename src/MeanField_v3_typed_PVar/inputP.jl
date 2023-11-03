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
num_laps = 2
num_ions = 4
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
tso = T0/40
tPh = 0.0*T0*ones(2)
tsoPh = 0.5*tso
model = FranssonHamiltonian
modelStr = "Fransson"

# Simulation Parameters
dt = 0.0002
tend = 15
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
#    phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)
popPhIn = zeros(N)
#xnIn,pnIn = [0.76,-3.61,-5.24,2.4],[-7.40,-4.33,1.59,-0.3]
xnIn,pnIn = zeros(N), zeros(N)
## Leads
P = 0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)

V = 60 * 0.2 * T0
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
path = "L$(num_laps)x$(num_ions)_tsoPh_sign"
if !isdir(path)
	mkpath(path)
end


num_avg = 1
for ka in 1:num_avg
    
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_2_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
    end
     
     
    dataname2 = "$(dataname1)_$(ka)"

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    SaveParams(params,parampath,purpose)
    
    if phIn == "Wigner"
	    	HamPar.X  = readdlm("init_X_traj$(version).dat")[ka,:]
		HamPar.P  = readdlm("init_X_traj$(version).dat")[ka,:]
    end
#    HamPar.X = xnIn
#    HamPar.P = pnIn
     
    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,elIntegrator,clIntegrator,popElIn)
    Average(ka,dataname1,num_avg)
end
