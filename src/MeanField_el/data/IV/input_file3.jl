module InputFile

#src_path=expanduser("~/Libraries/CISS/src/MeanField/")
#include("$(src_path)hamiltonian.jl")
#include("$(src_path)iosave.jl")
#include("$(src_path)multiopenint.jl")
#include("$(src_path)classical_phonons.jl")
include("../../hamiltonian.jl")
include("../../iosave.jl")
include("../../multiopenint.jl")

using .Hamiltonian
using .IOsave
using LinearAlgebra
using .MultiOpenIntegrator
# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 1500.0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

# Geometry Parameters
t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/4
num_laps = 3
ions = 4
l = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(l)*1.4*t0
qs = zeros(l)
H = Hamiltonian.FranssonHamiltonian(l,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)
if initProp
    #rhoIn = "pure"
    rhoIn = "equilibriumOpen"
    phIn = "TDSCF"
    # phIn = "Wigner"
    pathIn = ""
else
    rhoIn = "fromFile"
    phIn = rhoIn
    pathIn = "data/L3_dt0.001_1_final_Quantities"
end

popElIn = [1,1,0,0,0,0,0,0,0,0]

## Leads
tL = [t0/4,t0/4] 
tR = tL #[t/4,t/4] 
  #[1.2 for i in 1:numstep]
Vs = Array(0.1:0.1:2.5)
for V in Vs
	muL = V/2 #+ 6*t0
	muR = -V/2 #+ 6*t0
	nPoles = 25
	beta = 1/0.0001


	path = "V$(V)"
	if !isdir(path)
		mkdir(path)
	end

	num_avg = 1
	resampling = false
	rsmpl_stp = 15
    	dataname2 = "$(path)/L$(l)_laps$(num_laps)_ions$(ions)_V$(V)_tso$(tso)_beta$(beta)_dt$(abs(dt))"

	params = [("dt",dt),("Nstep",numstep),("N",l),("t0",t0),("tso",tso),("tph",tph),("tsoph",tsoph),("ws",ws),("tL",tL),("tR",tR),("V",V[1]),("nPoles",nPoles),("beta",beta)]

    	parampath = "$(path)/paramfile.txt"
    	purpose = "Test of Fransson model."
    	IOsave.SaveParams(params,parampath,purpose)

    	@timev Propagate(l,H,dt,numstep,tL,tR,beta,nPoles,muL,muR,dataname2,rhoIn,popElIn)
end

end

