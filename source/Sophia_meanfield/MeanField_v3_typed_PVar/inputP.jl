

include("hamiltonian.jl")
include("iosave.jl")
include("multiopenint.jl")
include("classical_phonons.jl")
include("structures.jl")
include("electronic_integration.jl")

using LinearAlgebra
using DelimitedFiles
# -------------- Parameters --------------------------


# Rudis Hamiltonian Parameters
num_laps = 1
num_ions =1
N =2* num_laps * num_ions # total number of sites in both chains --> * 2
ws  = 0.02*ones(N)
Hop = 0.1
T0  = Hop/ws[1]  # in eV
t0  = -T0
eu  = 0*T0
ed  = 3*T0
E   = [eu,ed]
gPh = zeros(2)
a   = 1.0
c   = 1.0
tso = 0.1*T0
tPh = 0.1*T0*ones(2)
tsoPh = 0.1*tso

#  veränderte und neue Variablen
#

onsight_e= E # Array([1. , 2.])./ ws[1] #j = 1; j = 2
thopping= Array([1.2*T0,-T0]) #j = 1, j= 2
lambda = -3*T0 #const.
t_so = tso # 1.2/ ws[1] #const
sites_per_turn = num_ions # man könnte noch überall die Variable sites_per... durch num_ions ersetzen
theta = pi/2

g_coupling =0.0# 0.05*T0  # = gPH ,  assumed constant
model = H_total_hampar
modelStr = "Guo"

# Simulation Parameters
dt = 0.0001/T0
tend = 3*dt#40/T0 # tend davor war 18 /T0 und hat 31 sec gebraucht
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


popElIn = zeros(2*N)
popElIn[1] = 1.0
popPhIn = zeros(N)
#xnIn,pnIn = [0.76,-3.61,-5.24,2.4],[-7.40,-4.33,1.59,-0.3]
#xnIn,pnIn = zeros(N), zeros(N)



xnIn,pnIn = ones(N), ones(N)

#xnIn, pnIn  = vcat( ones(Integer(N/2))*.5, ones(Integer(N/2))*1.5) ,ones(N)
## Leads
P =0.5
tl = 0#T0*10
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)
V= 0.
nPoles = 2
beta = ws[1]/0.025 
# hampar habe ich auch hier angepasst
HamPar = HamiltonianParams(
                        onsight_e,
                        thopping,
                        lambda,
                        t_so,
                        sites_per_turn,
                        theta,
                        g_coupling, # assumed const

                        num_laps,
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
path = "data$(version)M_test_g" # sd347@tcmp02:~/Bachelorarbeit/code/CISS/data
if !isdir(path)
        mkpath(path)
end

num_avg = 1
for ka in 1:num_avg # num_avg ist Anzahl propagationen für gleiche Zeit zum averagen

    if modelStr == "Fransson" # variablen aus hampar
        dataname1 = "$(path)/L$(N)_tPh_$(tPh[1])_thPd_$(tPh[2])_SOC$(tso)_SOCPh$(tsoPh)_V$(V)_w0$(ws[1])_P$(P)_$(clIntegrator)_$(elIntegrator)_2_dt$(abs(dt))"
    elseif modelStr == "Holstein"
        dataname1 = "$(path)/L$(N)_gPhU$(gPh[1])_gPhD$(gPh[2])_V$(V)_w0$(ws[1])_dt$(abs(dt))"
    elseif modelStr == "Guo"# hier Veränderungg
        dataname1 = "$(path)/L$(N)_g_coupling$(g_coupling)_V$(V)_w0$(ws[1])_num_laps$(num_laps)_num_ions$(num_ions)_dt$(abs(dt))" #???????????????????
    end


    dataname2 = "$(dataname1)_$(ka)"


        # true = propagate as usuals
        # false = use initial density matrix from previous computation

        initProp  = true # ???????????????? soll das dann auch für die erste Propagation auf true sein?
        if initProp
                Lauf = 0

                rhoIn = "pure"
             #   rhoIn = "equilibrium"
                phIn = "TDSCF"
        #	phIn = "Wigner"
                pathIn = ""
        else
                Lauf = 1

                rhoIn = "fromFile"
                phIn = rhoIn
                pathIn = "$(dataname2)_final_Quantities"
                 #"L1x2/L4_g_coupling0.39999999999999997_V11.999999999999998_w00.05_num_laps1_num_ions2_dt0.00125_1_final_Quantities"
        end

    params = [("dt",dt),("Nstep",numstep),("N",1),("ws",ws),("gPh",gPh),("tL",tL),("tR",tR),("V",V),("nPoles",nPoles),("beta",beta)]

    parampath = "$(path)/paramfile.txt"
    purpose = "Test of Guo model."
    SaveParams(params,parampath,purpose)

    if phIn == "Wigner"
                # das ist der Pfad zu daten aus generate_traj
                HamPar.X  = readdlm("init_X_traj$(version).dat")[ka,:]
                HamPar.P  = readdlm("init_P_traj$(version).dat")[ka,:]
        elseif phIn == "fromFile"
                file = h5open("$(pathIn)", "r")
                HamPar.X = read(file["X"])
                HamPar.P = read(file["P"])
        else
                #HamPar.X = xnIn
                #HamPar.P = pnIn

        end

    @timev Propagate(HamPar,dt,numstep,tL,tR,beta,nPoles,V,dataname2,rhoIn,elIntegrator,clIntegrator, initProp, Lauf, popElIn,pathIn)
    #Average(ka,dataname1,num_avg)
end

