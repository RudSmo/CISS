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
num_ions = 5
N = num_laps * num_ions
ws  = 0.2*0.002*ones(N)
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
dt = 0.00005
tend = 120
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

popElIn = zeros(2*N)
popPhIn = zeros(N)
#xnIn,pnIn = [0.7691068624482689,-3.61069243320506,-5.245657309206734,],[-7.401418656594664,-4.3318963287495675,1.5952234947490966]
xnIn,pnIn = zeros(N), zeros(N)
## Leads
P = -0.5
tl = T0/4
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)

v = 60
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
path = "L$(num_laps)x$(num_ions)"
if !isdir(path)
        mkpath(path)
end


num_avg = 10
for ka in 1:num_avg

    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
    end
    dataname2 = "$(dataname1)_$(ka)"

    if initProp
        rhoIn = "pure"
    #    rhoIn = "equilibrium"
    #    phIn = "TDSCF"
        phIn = "Wigner"
        pathIn = ""
        propStp = 0
    else
        rhoIn = "fromFile"
        phIn = rhoIn
        pathIn = "$(dataname2)_final_Quantities"
        propStp = 1
    end
    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Fransson model."
    if propStp == 0
            SaveParams(params,parampath,purpose)
    end

    if phIn == "Wigner"
                HamPar.X  = readdlm("init_X_traj$(version).dat")[ka,:] 
                HamPar.P  = readdlm("init_P_traj$(version).dat")[ka,:] 
		HamPar.X  = HamPar.X + 0.1 * rand(length(HamPar.X))
		HamPar.P  = HamPar.P + 0.1 * rand(length(HamPar.P))
    elseif phIn == "fromFile"
                file = h5open(pathIn)
                HamPar.X  = read(file["X"])
                HamPar.P  = read(file["P"])
    end

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,elIntegrator,clIntegrator,initProp,rhoIn,popElIn,pathIn,propStp)
    Average(ka,dataname1,num_avg)
end
