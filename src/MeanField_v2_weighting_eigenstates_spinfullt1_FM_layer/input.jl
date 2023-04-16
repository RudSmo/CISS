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
num_fm = 4
num_laps = 2
num_ions = 2
N = num_fm + num_laps * num_ions
fm_sites = Int.(1:num_fm)
mol_sites = Int.(num_fm+1:num_laps*num_ions)
coup_site = num_fm
T0  = 0.04
t0  = -T0
t0FM = -T0
coup_en = [-T0 0; 0 -T0]
eu  = 0.0 #-6*T0
ed  = 0.0 #-6*T0
E   = [eu,ed]
EFM = [eu,ed]
ws  = 0.01*ones(N)
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = T0/40
tPh = [-0.004,-0.004]
tsoPh = tso
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
#model = FranssonHamiltonianNH
sz = [0,0,1.]
mag = [sz for i in 1:num_fm]
Jsd = 1.0


# Simulation Parameters
dt = 0.002/T0
tend = 200/T0
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
    pathIn = "Holstein_L2/noEps/L2_gPhU0.004_gPhD0.004_V0.0_w00.002_X0_0.0_dt0.05_1_final_Quantities"
end
popElIn = zeros(2*N)
popPhIn = zeros(N)
vSampl = "v2"
phases = [rand()*2*pi for i in 1:length(popPhIn)]
x0 = 0.0
xnIn,pnIn = [(-1)^i*0.1*sin(0.001) for i in 1:N],[0.01*cos(0.001) for i in 1:N]


## Leads
P = 0.0
tl = 0.004
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)
V  = 0.0
nPoles = 10
beta = 1/0.025 
HamPar = HamiltonianParams( fm_sites,
			    mol_sites,
			    num_laps,
                            num_ions,
                            N,
                            model,
                            modelStr,
			    EFM,
                            E,
			    coup_site,
			    coup_en,
                            a,
                            c,
			    t0FM,
                            t0,
                            tPh,
                            tso,
                            tsoPh,
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
			    tL,
			    tR,
			    mag,
			    Jsd)



path = "L$(N)_nFM$(num_fm)_test"
if !isdir(path)
	mkpath(path)
end


num_avg = 1
for ka in 1:num_avg
    if modelStr == "Fransson"
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_noEps_dt$(abs(dt))"
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

