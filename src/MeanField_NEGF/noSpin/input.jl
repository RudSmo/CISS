include("hamiltonian.jl")
include("observables.jl")
include("structures.jl")
using PyPlot
using LinearAlgebra
using .Hamiltonian
using .Observables
using .CustomStructs

t0 = 1.0
eu = 0.0
ed = 0.0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/4
num_laps = 3
ions = 4
N = Int(num_laps * ions)
tph = 0.0
tsoph = 0.0
ws = ones(N)*1.4*t0
qs = zeros(N)

HS = Hamiltonian.FranssonHamiltonian(N,E,t0,tso,a,c,num_laps,tph,tsoph,ws,qs)


En = Array(-4.0:0.01:4.0)
sys_par = SystemParams(N,HS,En)

gLu = 0.25
gLd = 0.25
gLud = 0.0#1im*20.0*1e-1
gLdu = gLud

gRu = gLu
gRd = gLd
gRud = 0.0#1im*20.0*1e-1
gRdu = gRud

muL = 0.0
muR = 0.0
tL  = 0.0001
tR  = 0.0001

lead_par = LeadParamsSF(muL,muR,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tL,tR)


Vnu = Array(-5.0:0.03:5.0)
Vnd = Array(-5.0:0.03:5.0)


path = "data/Fransson_L$(N)"
mkpath(path)

filname = "$(path)/N$(N)_currents_spin.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vnu,Vnd,filname,savestp)
#end

