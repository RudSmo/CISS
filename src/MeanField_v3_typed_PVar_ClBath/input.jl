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
num_ions = 2
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
tPh = 0.1*T0*ones(2)
tsoPh = 0.0*tso
model = FranssonHamiltonian
modelStr = "Fransson"
ws = ws/2

# Auxiliary oscillators
maux = 50 # number of aux oscillators
eta = 0.05
wc = ws[1] # central frequency


w_rng = 1e-7:1e-5:0.015
ws_aux1 = w_rng[1]:(wc/2-w_rng[1])/(maux/4):wc/2
ws_aux2 = (wc/2):(3*wc/2)/(maux/2-1):2*wc
ws_aux3 = (2*wc):(w_rng[end]-3*wc/2)/(maux/4-1):w_rng[end]
ws_aux  = collect(Iterators.flatten([Array(ws_aux1),Array(ws_aux2),Array(ws_aux3)]))
gaux    = [sqrt(pi * eta * w * exp(-w/wc)/2) for w in ws_aux]  # Ohmic spectral density


# Simulation Parameters
dt = 0.0001
tend = 20
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
    rhoIn = "pure"
#    rhoIn = "equilibrium"
#    phIn = "TDSCF"
    phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)
popElIn[1] = 1.0
popPhIn = zeros(N)
xnIn,pnIn = zeros(N), zeros(N)
Xaux0,Paux0 = zeros(maux),zeros(maux)
## Leads
P = 0.0
tl = 0.0 * T0/4
tL = tl * [1,0]#(ones(2) .+ P * [1,-1])  
tR = tl * [1,0]#ones(2)

V = 75 * 0.2 * T0
nPoles = 1 #15
beta = ws[1]/0.025 

HamPar = HamiltonianParams( num_laps,
                            num_ions,
                            N,
			    maux,
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
			    gaux,
                            ws,
			    ws_aux,
                            xnIn,
			    Xaux0,
                            pnIn,
			    Paux0,
			    tL,
			    tR)


version = 0
path = "L$(num_laps)x$(num_ions)_Peierls_4MCTDH_comparison_V$(V)"
if !isdir(path)
	mkpath(path)
end


num_avg = 50
for ka in 1:num_avg
    
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_Naux$(maux)_2_dt$(abs(dt))"
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
		HamPar.Xaux = readdlm("init_Xaux_traj$(version).dat")[ka,:] 
		HamPar.Paux = readdlm("init_Paux_traj$(version).dat")[ka,:]
    end
#    HamPar.X = xnIn
#    HamPar.P = pnIn
     
    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,elIntegrator,clIntegrator,popElIn)
    Average(ka,dataname1,num_avg)
end
