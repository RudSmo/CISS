p = expanduser("~/Libraries/CISS/src/MeanField_NEGF/noSpin/") 
include("$(p)hamiltonian.jl")
include("$(p)observables.jl")
include("$(p)structures.jl")
using .Hamiltonian
using .Observables
using .CustomStructs

t0 = 0.04
eu = -6*t0
ed = -6*t0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/40
num_laps = 3
ions = 4
N = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(N)*1.4*t0
qs = zeros(N)

HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)


En = Array(-2.0:0.0005:1.0)
sys_par = SystemParams(N,HS,En)

P = 0.0
Gamma = t0/4
gLu = Gamma*(1 + P)
gLd = Gamma*(1 - P) 
gLud = 0.0
gLdu = gLud

gRu = Gamma
gRd = Gamma
gRud = 0.0
gRdu = gRud

Temp = 0.026
muL = 0.0
muR = 0.0
tL  = Temp
tR  = Temp

lead_par = LeadParamsSF(muL,muR,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tL,tR)


Vn = Array(0.0:0.01:2.5)


filname = "fS$(ions)_nL$(num_laps)_tso$(tso)_P$(P)_currents.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vn,filname,savestp)

