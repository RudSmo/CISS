module InputFile

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")
include("structures.jl")
using .Hamiltonian
using .IOsave
using LinearAlgebra
using DelimitedFiles
using .MultiOpenIntegrator
using .ClassicalPhonons
using .CustomStructures
# -------------- Parameters --------------------------


# Hamiltonian Parameters
num_laps = 1
num_ions = 3
N = num_laps * num_ions
T0  = 0.04
t0  = 0.0
eu  = 0.0
ed  = 0.0
E   = [eu,ed]
ws  = 0.002*ones(N)
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.0
tPh = [-0.004,-0.004]
tsoPh = 0.0001
#model = HolsteinHopping
#modelStr = "Holstein"
model = FranssonHamiltonian
modelStr = "Fransson"
#model = FranssonHamiltonianNH

# Simulation Parameters
dt = 0.002/T0
tend = 1000/T0
numstep = Int(tend/abs(dt))
writeMode = "minimal"
initProp  = true
# true = propagate as usuals
# false = use initial density matrix from previous computation

if initProp
    #rhoIn = "pure"
    rhoIn = "equilibrium"
    phIn = "TDSCF"
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
xnIn,pnIn = [1/sqrt(2*ws[1]) , 0],zeros(N)


## Leads
P = 0.0
tl = 0.004
tL = tl * (ones(2) .+ P * [1,-1])  
tR = tl * ones(2)
V  = 0.2 * T0 * 65
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
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
			    tL,
			    tR)

H = model(HamPar)
path = "L$(num_laps)x$(num_ions)_hamiltonian"
if !isdir(path)
       mkpath(path)
end
fil_R = open("$(path)/hamiltonian_$(num_laps)x$(num_ions)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_EL_real.dat", "w")
fil_I = open("$(path)/hamiltonian_$(num_laps)x$(num_ions)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_EL_imag.dat", "w")
writedlm(fil_R,real.(H))
writedlm(fil_I,imag.(H))


#path = "L$(num_laps)x$(num_ions)_spinSoliton"
#if !isdir(path)
#	mkpath(path)
#end
#
#
#num_avg = 1
#for ka in 1:num_avg
#    if modelStr == "Fransson"
#        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_noEps_dt$(abs(dt))"
#    elseif modelStr == "Holstein"
#        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_X0_$(x0)_dt$(abs(dt))"
#    end
#     
#     
#    dataname2 = "$(dataname1)_$(ka)"
#
#    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]
#
#    parampath = "$(path)/paramfile.txt"
#    purpose = "Test of Fransson model."
#    IOsave.SaveParams(params,parampath,purpose)
#
#    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,phIn,vSampl,popElIn)
##    Average(ka,dataname1,num_avg)
#end

end

